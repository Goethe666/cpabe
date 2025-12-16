#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from charm.toolbox.pairinggroup import PairingGroup, GT, G1, G2, ZR
from MSP import MSP
from revocation_tree import CompleteBinaryTree as RevCompleteBinaryTree
import hashlib


class IndirectRevocationABE:
    """基于二叉树的纯间接撤销ABE系统"""
    
    def __init__(self, num_users=16, group_type='MNT224'):
        """
        初始化间接撤销ABE系统
        
        Args:
            num_users (int): 最大用户数量
            group_type (str): 配对群类型
        """
        
        # 初始化密码学组件
        self.group = PairingGroup(group_type)
        self.msp = MSP(self.group, verbose=False)
        
        # 初始化二叉树（引用外部撤销树实现）
        self.tree = RevCompleteBinaryTree(num_users)
        self.num_users = num_users
        
        # 系统参数（不在对象中持有密钥）
        self.current_time = 1
        
        # 存储（不保存私钥与策略，按需传参）
        self.revoked_users = set()  # 撤销用户集合
        
        # 不在构造函数中执行设置；密钥通过 setup() 返回

    def setup(self):
        return self._setup()
        
    
    def _setup(self):
        """系统设置算法"""
        
        # 生成主密钥 α
        alpha = self.group.random(ZR)
        
        # 生成生成元（非对称）：G1 与 G2
        g = self.group.random(G1)
        g2 = self.group.random(G2)
        
        # 生成随机元素 u_0, ..., u_d, h_0, ..., h_m
        # 这里简化为固定数量
        u_elements = [self.group.random(G1) for _ in range(10)]
        h_elements = [self.group.random(G1) for _ in range(10)]
        
        # 设置公钥
        # 预计算 GT 生成元 e(g1,g2)
        e_gg = self.group.pair_prod(g, g2)
        public_key = {
            'g': g,
            'g2': g2,
            'e_gg': e_gg,
            'e_g_g_alpha': e_gg ** alpha,
            'u_elements': u_elements,
            'h_elements': h_elements
        }
        
        # 为树中每个节点生成多项式系数
        node_coeffs = {}
        for node_id, node in self.tree.nodes.items():
            a_x = self.group.random(ZR)
            node.polynomial_coeff = a_x
            node_coeffs[node_id] = a_x
        
        # 设置主密钥
        master_key = {
            'alpha': alpha,
            'node_coeffs': node_coeffs
        }

        return public_key, master_key
        
    
    def _hash_to_zp(self, value_str, value_type):
        """
        [辅助函数] 将字符串安全地哈希到 Zp (指数域)
        
        Args:
            value_str (str): 要哈希的原始字符串
            value_type (str): 值的类型
            
        Returns:
            ZR: 在 Zp 中的一个元素
        """
        
        prefix = ''
        
        # 1. 根据类型选择前缀，严格遵循论文 的建议
        if value_type == 'attribute':
            # 对应论文中的 H(0||x)
            prefix = '0' 
        elif value_type == 'time':
            # 对应论文中的 H(1||t)
            prefix = '1'
        else:
            # 如果类型错误，立即报错
            raise ValueError(f"未知的哈希类型: {value_type}")

        # 2. 将前缀和值连接起来
        data = (prefix + str(value_str)).encode()
        
        # 3. 使用完整的 SHA-256 哈希 (不截断)
        full_hash = hashlib.sha256(data).hexdigest()
        
        # 4. 将完整哈希转换为大整数
        hash_int = int(full_hash, 16)
        
        # 5. 对群的阶 p (即 group.order()) 取模，并返回一个 Zp 元素
        return self.group.init(ZR, hash_int % self.group.order())


    def _get_poly_value(self, master_key, node, z_value):
        """
        计算多项式 f_x(z) = a_x * z + α 的值
        
        Args:
            node: 树节点，包含 polynomial_coeff (a_x)
            z_value: 多项式的输入值 z (可以是整数或ZR元素)
            
        Returns:
            ZR: f_x(z) 的计算结果
        """
        # 确保 z_value 是 ZR 元素
        if isinstance(z_value, int):
            z_in_zp = self.group.init(ZR, z_value)
        else:
            z_in_zp = z_value
            
        # 计算 f_x(z) = a_x * z + α
        return node.polynomial_coeff * z_in_zp + master_key['alpha']

    def function_F(self, public_key, x):    
        """
        哈希函数F：实现论文中的 F(x) = prod(h_j ^ (x^j))
        (F函数用于普通属性 ω)
        """
        
        # 步骤 1: 将属性 "ω" (即 x) 哈希到 Zp
        # (x_in_zp 就是公式 中的 x)
        x_in_zp = self._hash_to_zp(x, value_type='attribute')    
        
        h = public_key['h_elements'] 
        m = len(h) - 1
        
        # 在指数上计算多项式
        F_x = h[0]
        x_pow = x_in_zp
        
        for j in range(1, m + 1):
            term = h[j] ** x_pow
            F_x = F_x * term
            x_pow = x_pow * x_in_zp
            
        return F_x

    def function_P(self, public_key, x):
        """
        哈希函数P：实现论文 中的 P(x) = prod(u_j ^ (x^j))
        (P函数在间接撤销中用于 "时间 t")
        
        Args:
            x (str or int): 要哈希的时间 (例如 1, "2025-Q1")
        """
        
        # 步骤 1: 将时间 "t" (即 x) 哈希到 Zp
        # (x_in_zp 就是公式 中的 x)
        x_in_zp = self._hash_to_zp(x, value_type='time')    
        
        u = public_key['u_elements'] # u_0, u_1, ...
        d = len(u) - 1
        
        # 步骤 2: 在指数上计算多项式
        
        # 对应 j=0, 即 u_0
        P_x = u[0] 
        x_pow = x_in_zp # 准备 j=1 时的 x^1
        
        # 对应 j=1 到 d
        for j in range(1, d + 1):
            term = u[j] ** x_pow  # u_j^(x^j)
            P_x = P_x * term      # ∏ (连乘)
            
            # 为下一次循环 j+1 准备 x^(j+1)
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

            for node in path:
                f_x_1 = self._get_poly_value(master_key, node, self.group.init(ZR, 1))

                v = [f_x_1]
                for _ in range(num_cols - 1):
                    v.append(self.group.random(ZR))

                node_key_shares = {}
                for attr_with_index, row in mono_msp.items():
                    share_lambda = self.group.init(ZR, 0)
                    cols_to_mult = len(row)
                    for i in range(cols_to_mult):
                        share_lambda += row[i] * v[i]

                    attr_plain = attr_map[attr_with_index]
                    r_xi = self.group.random(ZR)
                    F_attr = F_cache[attr_plain]
                    D_1 = (public_key['g'] ** share_lambda) * (F_attr ** r_xi)
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
        """
        加密算法（修复代数错误版本）
        
        Args:
            message_element (GT): 已映射到GT群的消息元素
            attributes (list): 属性列表
            time_period (int): 时间周期
            
        Returns:
            dict: 密文字典，包含C, C1, C2, C3等组件
        """
        try:
            
            # 生成随机数 s
            s = self.group.random(ZR)
            
            # 计算 C = m * (e(g,g)^α)^s
            # 修复：直接使用 message_element，避免 ZR * GT 的类型错误
            C = message_element * (public_key['e_g_g_alpha'] ** s)
            
            # 计算 C^{(1)} = g2^s（置于 G2 以适配非对称配对）
            C1 = public_key['g2'] ** s
            
            # 计算 C_k^{(2)} = F(k)^s 对每个属性 k
            C2 = {}
            for attr in attributes:
                F_attr = self.function_F(public_key, attr)
                C2[attr] = F_attr ** s
            
            # 计算 C^{(3)} = P(t)^s
            P_t = self.function_P(public_key, time_period)
            C3 = P_t ** s
            
            # 构造密文
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
        """
        生成密钥更新材料（支持间接撤销机制）
        
        Args:
            time_period (int): 时间周期
            
        Returns:
            dict: 密钥更新材料，包含覆盖节点的更新组件
        """
        try:
            
            # 使用覆盖算法找到覆盖所有未撤销用户的最小节点集
            cover_nodes = self.tree.cover_algorithm(self.revoked_users)
            
            # 为每个覆盖节点生成更新材料
            update_material = {}
            
            # 统一将时间哈希到 Zp
            t_zr = self._hash_to_zp(time_period, 'time')
            
            for node in cover_nodes:
                # 修复：重新引入随机化 r_prime 以确保解密正确性
                r_prime = self.group.random(ZR)
                
                # 计算 f_x(t) 使用哈希到 ZR 的时间值
                f_x_t = self._get_poly_value(master_key, node, t_zr)
                
                # 计算 U_x^{(1)} = g^{f_x(t)} * P(t)^{r'}
                P_t = self.function_P(public_key, time_period)
                U_1 = (public_key['g'] ** f_x_t) * (P_t ** r_prime)
                
                # 计算 U_x^{(2)} = g2^{r'}（置于 G2 以适配非对称配对）
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
        """
        核心解密算法（修复LSSS逻辑问题版本）
        
        Args:
            ciphertext (dict): 密文
            user_key (dict): 用户私钥
            update_material (dict): 密钥更新材料
            user_id (int): 用户ID
            
        Returns:
            GT: 解密后的消息元素，失败返回None
        """
        try:
            
            # 检查用户是否被撤销
            if user_id in self.revoked_users:
                return None
            
            # 获取用户路径
            path = self.tree.get_path(user_id)
            if not path:
                return None
            
            # 找到用户路径上在更新材料中的节点
            valid_node = None
            for node in path:
                if node.node_id in update_material:
                    valid_node = node
                    break
            
            if valid_node is None:
                return None
            
            
            # 获取节点密钥和更新材料
            node_user_key = user_key[valid_node.node_id]
            node_update = update_material[valid_node.node_id]

        
            # 1. 从用户密钥中获取访问策略（在 keygen 时已存储）
            user_policy_obj = node_user_key.get('policy')
            user_policy_str = node_user_key.get('policy_str')
            if not user_policy_obj or not user_policy_str:
                return None

            # 2. 从密文中获取属性列表
            ciphertext_attributes = ciphertext.get('attributes')
            if not ciphertext_attributes:
                return None

            # print(f"  用户策略: {user_policy_str}")
            # print(f"  密文属性: {ciphertext_attributes}")

            # 3. LSSS 步骤：剪枝 + 系数计算（直接使用 Charm 底层接口，避免类型混淆）
            pruned_list = self.msp.prune(user_policy_obj, ciphertext_attributes)
            if not pruned_list:
                return None

            coeffs_all = self.msp.getCoefficients(user_policy_obj)

            # 4. 将系数保留为带索引的键，匹配 keygen 存储的主键
            coeffs_indexed = {}
            attrs_plain_for_log = []
            for node in pruned_list:
                attr_idx = node.getAttributeAndIndex()
                coeff = coeffs_all.get(attr_idx, None)
                if coeff is None:
                    continue
                coeffs_indexed[attr_idx] = coeff
                try:
                    attrs_plain_for_log.append(self.msp.strip_index(attr_idx))
                except Exception:
                    attrs_plain_for_log.append(attr_idx)

            # print(f"  满足策略的属性(去索引显示): {attrs_plain_for_log}")
            
            # 计算 K' = ∏ (e(D_{x,i}^{(1)}, C^{(1)}) / e(C_k^{(2)}, D_{x,i}^{(2)}))^{ν_i}
            K_prime = self.group.init(GT, 1)  # 初始化为GT群的单位元
            
            for attr_idx, coeff in coeffs_indexed.items():
                if attr_idx in node_user_key['node_key']:
                    # 获取密钥组件（带索引键）
                    D_1 = node_user_key['node_key'][attr_idx]['D_1']
                    D_2 = node_user_key['node_key'][attr_idx]['D_2']
                    # 获取用于查找 C2 的纯属性名
                    attr_plain = node_user_key['node_key'][attr_idx].get('attr_plain', attr_idx)
                    if attr_plain in ciphertext['C2']:
                        # 获取密文组件
                        C1 = ciphertext['C1']
                        C2_attr = ciphertext['C2'][attr_plain]
                        # 计算配对
                        numerator = self.group.pair_prod(D_1, C1)
                        denominator = self.group.pair_prod(C2_attr, D_2)
                        # 计算 (e(D_1, C1) / e(C2, D_2))^coeff
                        term = (numerator / denominator) ** coeff
                        K_prime = K_prime * term
            
            # 计算 K_{time} = e(U^{(1)}, C^{(1)}) / e(C^{(3)}, U^{(2)})
            U_1 = node_update['U_1']
            U_2 = node_update['U_2']
            C1 = ciphertext['C1']
            C3 = ciphertext['C3']
            
            numerator_time = self.group.pair_prod(U_1, C1)
            denominator_time = self.group.pair_prod(C3, U_2)
            K_time = numerator_time / denominator_time
            
            # 拉格朗日插值：K = (K')^{t/(t-1)} * (K_{time})^{1/(1-t)}
            t = ciphertext['time_period']
            # 统一使用哈希后的时间值进行计算
            t_zr = self._hash_to_zp(t, 'time')
            
            # 计算 t/(t-1)
            # 确保 t 的哈希值不为 1
            if t_zr == 1:
                return None
            
            lagrange_coeff1 = t_zr / (t_zr - 1)
            
            # 计算 1/(1-t)
            lagrange_coeff2 = 1 / (1 - t_zr)
            
            # 计算最终密钥 K
            K = (K_prime ** lagrange_coeff1) * (K_time ** lagrange_coeff2)
            
            # 恢复消息：m = C / K
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
