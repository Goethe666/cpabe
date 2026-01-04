
from charm.toolbox.pairinggroup import PairingGroup, GT, G1, G2, ZR
from MSP import MSP
from revocation_tree import CompleteBinaryTree as RevCompleteBinaryTree
import hashlib

class IndirectRevocationABE:
    def __init__(self, num_users=16, group_type='MNT224'):
        self.group = PairingGroup(group_type)
        self.msp = MSP(self.group, verbose=False)
        self.tree = RevCompleteBinaryTree(num_users)
        self.num_users = num_users
        self.current_time = 1
        self.revoked_users = set()
        self._hash_cache = {}

    def setup(self):
        return self._setup()
        
    
    def _setup(self):
        alpha = self.group.random(ZR)
        g = self.group.random(G1)
        g2 = self.group.random(G2)
        u_elements = [self.group.random(G1) for _ in range(10)]
        h_elements = [self.group.random(G1) for _ in range(10)]
        e_gg = self.group.pair_prod(g, g2)
        public_key = {
            'g': g,
            'g2': g2,
            'e_gg': e_gg,
            'e_g_g_alpha': e_gg ** alpha,
            'u_elements': u_elements,
            'h_elements': h_elements
        }
        
        node_coeffs = {}
        for node_id, node in self.tree.nodes.items():
            a_x = self.group.random(ZR)
            node.polynomial_coeff = a_x
            node_coeffs[node_id] = a_x
        master_key = {
            'alpha': alpha,
            'node_coeffs': node_coeffs
        }

        return public_key, master_key
        
    
    def _hash_to_zp(self, value_str, value_type):
        key = (value_type, str(value_str))
        if key in self._hash_cache:
            return self._hash_cache[key]
        prefix = ''
        if value_type == 'attribute':
            prefix = '0' 
        elif value_type == 'time':
            prefix = '1'
        else:
            raise ValueError(f"未知的哈希类型: {value_type}")
        data = (prefix + str(value_str)).encode()
        full_hash = hashlib.sha256(data).hexdigest()
        hash_int = int(full_hash, 16)
        res = self.group.init(ZR, hash_int % self.group.order())
        self._hash_cache[key] = res
        return res


    def _get_poly_value(self, master_key, node, z_value):
        if isinstance(z_value, int):
            z_in_zp = self.group.init(ZR, z_value)
        else:
            z_in_zp = z_value
        return node.polynomial_coeff * z_in_zp + master_key['alpha']

    def function_F(self, public_key, x):    
        x_in_zp = self._hash_to_zp(x, value_type='attribute')    
        h = public_key['h_elements'] 
        m = len(h) - 1
        F_x = h[0]
        x_pow = x_in_zp
        for j in range(1, m + 1):
            term = h[j] ** x_pow
            F_x = F_x * term
            x_pow = x_pow * x_in_zp
        return F_x

    def function_P(self, public_key, x):
        x_in_zp = self._hash_to_zp(x, value_type='time')    
        u = public_key['u_elements']
        d = len(u) - 1
        P_x = u[0] 
        x_pow = x_in_zp
        for j in range(1, d + 1):
            term = u[j] ** x_pow
            P_x = P_x * term
            x_pow = x_pow * x_in_zp
        return P_x
    def keygen(self, public_key, master_key, user_id, access_policy):
        try:
            if user_id not in self.tree.leaves:
                return False

            policy_obj = self.msp.createPolicy(access_policy)
            mono_msp = self.msp.convert_policy_to_msp(policy_obj)
            num_cols = self.msp.len_longest_row

            attr_map = {}
            F_cache = {}
            for attr in mono_msp.keys():
                plain_attr = self.msp.strip_index(attr)
                attr_map[attr] = plain_attr
                if plain_attr not in F_cache:
                    F_cache[plain_attr] = self.function_F(public_key, plain_attr)

            path = self.tree.get_path(user_id)
            user_key = {}

            mono_items = list(mono_msp.items())
            row_specs = {}
            for attr_with_index, row in mono_items:
                pos = []
                neg = []
                other = []
                for i, val in enumerate(row):
                    if val == 0:
                        continue
                    if val == 1:
                        pos.append(i)
                        continue
                    if val == -1:
                        neg.append(i)
                        continue
                    other.append((i, val))
                row_specs[attr_with_index] = (pos, neg, other)

            if not hasattr(self, '_pk_pp_done'):
                self._pk_pp_done = set()
            pk_id = id(public_key)
            if pk_id not in self._pk_pp_done:
                public_key['g'].initPP()
                public_key['g2'].initPP()
                self._pk_pp_done.add(pk_id)

            for node in path:
                f_x_1 = self._get_poly_value(master_key, node, 1)

                v = [f_x_1]
                for _ in range(num_cols - 1):
                    v.append(self.group.random(ZR))

                G_v = [public_key['g'] ** val for val in v]

                neg_union = set()
                for attr_with_index, _row in mono_items:
                    _pos, neg, _other = row_specs[attr_with_index]
                    for i in neg:
                        neg_union.add(i)
                inv_G_v = {}
                for i in neg_union:
                    inv_G_v[i] = G_v[i] ** -1

                node_key_shares = {}
                for attr_with_index, _row in mono_items:
                    D_1_base = self.group.init(G1, 1)
                    pos, neg, other = row_specs[attr_with_index]
                    for i in pos:
                        D_1_base *= G_v[i]
                    for i in neg:
                        D_1_base *= inv_G_v[i]
                    for i, val in other:
                        D_1_base *= (G_v[i] ** val)

                    attr_plain = attr_map[attr_with_index]
                    r_xi = self.group.random(ZR)
                    F_attr = F_cache[attr_plain]

                    D_1 = D_1_base * (F_attr ** r_xi)
                    D_2 = public_key['g2'] ** r_xi
                    node_key_shares[attr_with_index] = {
                        'D_1': D_1,
                        'D_2': D_2,
                        'attr_plain': attr_plain
                    }

                user_key[node.node_id] = {
                    'node_key': node_key_shares,
                    'policy': policy_obj,
                    'policy_str': access_policy
                }

            return user_key
        except Exception as e:
            return False

    def encrypt(self, public_key, message_element, attributes, time_period):
        try:
            s = self.group.random(ZR)
            C = message_element * (public_key['e_g_g_alpha'] ** s)
            C1 = public_key['g2'] ** s
            C2 = {}
            
            for attr in attributes:
                F_attr = self.function_F(public_key, attr)
                C2[attr] = F_attr ** s
                
            P_t = self.function_P(public_key, time_period)
            C3 = P_t ** s
            ciphertext = {
                'C': C,
                'C1': C1,
                'C2': C2,
                'C3': C3,
                'attributes': attributes,
                'time_period': time_period
            }
            
            return ciphertext
            
        except Exception as e:
             return None


    def key_update(self, public_key, master_key, time_period):
        try:
            cover_nodes = self.tree.cover_algorithm(self.revoked_users)
            update_material = {}
            t_zr = self._hash_to_zp(time_period, 'time')
            P_t = self.function_P(public_key, time_period)
            for node in cover_nodes:
                r_prime = self.group.random(ZR)
                f_x_t = self._get_poly_value(master_key, node, t_zr)
                U_1 = (public_key['g'] ** f_x_t) * (P_t ** r_prime)
                U_2 = public_key['g2'] ** r_prime
                
                update_material[node.node_id] = {
                    'U_1': U_1,
                    'U_2': U_2,
                    'time_period': time_period
                }
            
            return update_material
            
        except Exception as e:
            return None

    def _core_decrypt(self, ciphertext, user_key, update_material, user_id):
        try:
            
            if user_id in self.revoked_users:
                return None
            path = self.tree.get_path(user_id)
            if not path:
                return None
            valid_node = None
            for node in path:
                if node.node_id in update_material:
                    valid_node = node
                    break
            
            if valid_node is None:
                return None
            node_user_key = user_key[valid_node.node_id]
            node_update = update_material[valid_node.node_id]

        
            user_policy_obj = node_user_key.get('policy')
            if not user_policy_obj:
                return None
            ciphertext_attributes = ciphertext.get('attributes')
            if not ciphertext_attributes:
                return None
            pruned_list = self.msp.prune(user_policy_obj, ciphertext_attributes)
            if not pruned_list:
                return None
            coeffs_all = self.msp.getCoefficients(user_policy_obj)
            coeffs_indexed = {}
            for node in pruned_list:
                attr_idx = node.getAttributeAndIndex()
                coeff = coeffs_all.get(attr_idx, None)
                if coeff is None:
                    continue
                coeffs_indexed[attr_idx] = coeff
            K_prime = self.group.init(GT, 1)
            
            for attr_idx, coeff in coeffs_indexed.items():
                if attr_idx in node_user_key['node_key']:
                    D_1 = node_user_key['node_key'][attr_idx]['D_1']
                    D_2 = node_user_key['node_key'][attr_idx]['D_2']
                    attr_plain = node_user_key['node_key'][attr_idx].get('attr_plain', attr_idx)
                    if attr_plain in ciphertext['C2']:
                        C1 = ciphertext['C1']
                        C2_attr = ciphertext['C2'][attr_plain]
                        numerator = self.group.pair_prod(D_1, C1)
                        denominator = self.group.pair_prod(C2_attr, D_2)
                        term = (numerator / denominator) ** coeff
                        K_prime = K_prime * term
            U_1 = node_update['U_1']
            U_2 = node_update['U_2']
            C1 = ciphertext['C1']
            C3 = ciphertext['C3']
            numerator_time = self.group.pair_prod(U_1, C1)
            denominator_time = self.group.pair_prod(C3, U_2)
            K_time = numerator_time / denominator_time
            t = ciphertext['time_period']
            t_zr = self._hash_to_zp(t, 'time')
            if t_zr == 1:
                return None
            lagrange_coeff1 = t_zr / (t_zr - 1)
            lagrange_coeff2 = 1 / (1 - t_zr)
            K = (K_prime ** lagrange_coeff1) * (K_time ** lagrange_coeff2)
            message_element = ciphertext['C'] / K
            
            return message_element
            
        except Exception as e:
             return None

    
     
    
    def revoke_user(self, user_id, reason=""):
        """
        撤销用户
        
        Args:
            user_id (int): 用户ID
            reason (str): 撤销原因
        """
        if user_id in self.revoked_users:
            return
        
        self.revoked_users.add(user_id)
    
    def advance_time(self):
        """推进时间"""
        self.current_time += 1
        return self.current_time

    
    def demo_indirect_revocation(self):
        pk, mk = self.setup()

        sk1 = self.keygen(pk, mk, 1, "(EMPLOYEE and FINANCE)")
        sk2 = self.keygen(pk, mk, 2, "(EMPLOYEE and HR)")
        sk3 = self.keygen(pk, mk, 3, "(MANAGER and (FINANCE or HR))")
        sk4 = self.keygen(pk, mk, 4, "ADMIN")

        message1 = self.group.random(GT)
        message2 = self.group.random(GT)
        message3 = self.group.random(GT)

        ct1_data = self.encrypt(pk, message1, ["EMPLOYEE", "FINANCE"], 2)
        ct2_data = self.encrypt(pk, message2, ["EMPLOYEE", "HR"], 2)
        ct3_data = self.encrypt(pk, message3, ["MANAGER", "FINANCE"], 2)

        uk1 = self.key_update(pk, mk, 2)

        decrypted1 = self._core_decrypt(ct1_data, sk1, uk1, 1)
        decrypted2 = self._core_decrypt(ct2_data, sk2, uk1, 2)
        decrypted3 = self._core_decrypt(ct3_data, sk3, uk1, 3)

        print(f"用户1解密与明文相同: {decrypted1 is not None and decrypted1 == message1}")
        print(f"用户2解密与明文相同: {decrypted2 is not None and decrypted2 == message2}")
        print(f"用户3解密与明文相同: {decrypted3 is not None and decrypted3 == message3}")

        self.revoke_user(1, "违反公司政策")

        new_time = self.advance_time()
        uk2 = self.key_update(pk, mk, new_time)

        decrypted1_after = self._core_decrypt(ct1_data, sk1, uk2, 1)
        print(f"用户1撤销后解密与明文相同: {decrypted1_after is not None and decrypted1_after == message1}")

        message4 = self.group.random(GT)
        ct4_data = self.encrypt(pk, message4, ["MANAGER", "FINANCE"], new_time)

        decrypted4_user1 = self._core_decrypt(ct4_data, sk1, uk2, 1)
        decrypted4_user3 = self._core_decrypt(ct4_data, sk3, uk2, 3)
        print(f"用户1新消息解密与明文相同: {decrypted4_user1 is not None and decrypted4_user1 == message4}")
        print(f"用户3新消息解密与明文相同: {decrypted4_user3 is not None and decrypted4_user3 == message4}")

def main():
    """主函数"""
    try:
        # 创建间接撤销ABE系统
        abe_system = IndirectRevocationABE(num_users=8)
        
        # 运行演示
        abe_system.demo_indirect_revocation()
        
    except Exception as e:
        pass

if __name__ == "__main__":
    main()
