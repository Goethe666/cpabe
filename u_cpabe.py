from charm.toolbox.pairinggroup import PairingGroup, G1, G2, GT, ZR
from MSP import MSP
from revocation_tree import CompleteBinaryTree
import math

class UnifiedCPABE:
    def __init__(self, d=2, group_type='MNT224'):
        self.group = PairingGroup(group_type)
        self.d = d
        self.N = self.d * self.d
        self.tree = CompleteBinaryTree(self.N)
        self.node_coeffs = {}
        for node_id, node in self.tree.nodes.items():
            a_n = self.group.random(ZR)
            node.polynomial_coeff = a_n
            self.node_coeffs[node_id] = a_n
        self.msp = MSP(self.group, verbose=False)
        self.bHash = self.group.hash(str(int(self.group.order()) + 1), G1)
        self.revoked_users = set()

    def setup(self):
        d=self.d
        g1 = self.group.random(G1)
        g2 = self.group.random(G2)
        e_gg = self.group.pair_prod(g1, g2)
        alpha_rows = [self.group.random(ZR) for _ in range(d)]
        r_rows = [self.group.random(ZR) for _ in range(d)]
        c_cols = [self.group.random(ZR) for _ in range(d)]
        alpha = self.group.random(ZR)
        mpk = {
            'g1': g1,
            'g2': g2,
            'E': [g2 ** r for r in r_rows],
            'F': [g1 ** c for c in c_cols],
            'G': [e_gg ** a for a in alpha_rows],
            'e_gg_alpha': e_gg ** alpha
        }
        msk = {
            'alpha_rows': alpha_rows,
            'r_rows': r_rows,
            'c_cols': c_cols,
            'alpha': alpha,
            'node_coeffs': self.node_coeffs
        }
        return mpk, msk

    def _user_id(self, i, j):
        return (i - 1) * self.d + j

    def _H1(self, val):
        return self.group.hash(str(val), G1)

    def _H2(self, mpk, t):
        try:
            t_str = str(int(t))
        except Exception:
            t_str = str(t)
        return  mpk['g1']**self.group.hash(t_str, ZR)

    def _index_to_pos(self, k):
        if k < 1 or k > self.N:
            return None
        i = ((k - 1) // self.d) + 1
        j = ((k - 1) % self.d) + 1
        return (i, j)

    def keygen(self, mpk, msk, S, i, j):
        sigma = self.group.random(ZR)
        sk1 = (mpk['g1'] ** msk['alpha_rows'][i - 1]) * (mpk['g1'] ** (msk['r_rows'][i - 1] * msk['c_cols'][j - 1])) * (mpk['g1'] ** sigma)
        sk2 = mpk['g2'] ** sigma
        sk3 = {}
        uid = self._user_id(i, j)
        path = self.tree.get_path(uid)
        sk4 = {}
        sk5={}
        hashed_attrs = {u: self._H1(u) for u in S}
        for node in path:
            sigmax=self.group.random(ZR)
            sk3[node.node_id] = {u: hashed_attrs[u] ** (sigma+sigmax) for u in S}
            f1 = node.polynomial_coeff * self.group.init(ZR, 1) + msk['alpha']
            sk4[node.node_id] = (mpk['g1'] ** f1) * (self.bHash ** (sigma+sigmax))
            sk5[node.node_id] = (mpk['g2'] ** (sigma+sigmax))
        return {
            'sk1': sk1,
            'sk2': sk2,
            'sk3': sk3,
            'sk4': sk4,
            'sk5':sk5,
            'S': S,
            'pos': (i, j)
        }

    def encrypt(self, mpk, M, policy_str, t):
        policy = self.msp.createPolicy(policy_str)
        mono = self.msp.convert_policy_to_msp(policy)
        num_cols = self.msp.len_longest_row
        tau = self.group.random(ZR)
        eta = self.group.random(ZR)
        a = self.group.random(ZR)
        b = self.group.random(ZR)
        c = self.group.random(ZR)
        s_rows = [self.group.random(ZR) for _ in range(self.d)]
        w = {y: [self.group.random(ZR), self.group.random(ZR), self.group.random(ZR)] for y in range(1, self.d + 1)}
        v = [self.group.random(ZR) for _ in range(num_cols)]
        s = self.group.random(ZR)
        r=[self.group.random(ZR) for _ in range(self.d)]
        v[0] = s
        v1 = [a, self.group.init(ZR, 0), c]
        v2 = [self.group.init(ZR, 0), b, c]
        v3 = [v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0]]
        v_c = [self.group.random(ZR), self.group.random(ZR), self.group.random(ZR)]
        R = {}
        Rt = {}
        A = {}
        Apr = {}
        B = {}
        for x in range(1, self.d + 1):
            qx = self.group.random(ZR)
            qx2 = self.group.random(ZR)
            vx = [qx * v1[0] + qx2 * v2[0], qx * v1[1] + qx2 * v2[1], qx * v1[2] + qx2 * v2[2]]
            dot_v = vx[0] * v_c[0] + vx[1] * v_c[1] + vx[2] * v_c[2]
            R[x] = [mpk['E'][x - 1] ** (s_rows[x - 1] * vx[0]), mpk['E'][x - 1] ** (s_rows[x - 1] * vx[1]), mpk['E'][x - 1] ** (s_rows[x - 1] * vx[2])]
            Rt[x] = [R[x][0] ** eta, R[x][1] ** eta, R[x][2] ** eta]
            A[x] = mpk['g2'] ** (s_rows[x - 1] * tau * dot_v)
            Apr[x] =(mpk['g1'] ** (s_rows[x - 1] * tau * dot_v))
            B[x] = M * (mpk['G'][x - 1] ** (s_rows[x - 1] * tau * dot_v)) * (mpk['e_gg_alpha'] ** s)
        C = {}
        Ct = {}
        for y in range(1, self.d + 1):
            C[y] = [mpk['F'][y - 1] ** (tau * v_c[0]) * (mpk['g1']** (w[y][0] * eta)), mpk['F'][y - 1] ** ( tau * v_c[1]) * (mpk['g1'] ** (w[y][1] * eta)), mpk['F'][y - 1] ** ( tau * v_c[2]) * (mpk['g1'] ** (w[y][2] * eta))] 
            Ct[y] = [mpk['g1'] ** w[y][0], mpk['g1'] ** w[y][1], mpk['g1'] ** w[y][2]]
        Pl = {}
        Pr = {}
        s_prime = self.group.random(ZR)
        for attr, row in mono.items():
            attr_stripped = self.msp.strip_index(attr)
            attrHash = self.group.hash(attr_stripped, G1)
            len_row = len(row)
            Mivtop = self.group.init(ZR, 0)
            for i in range(len_row):
                Mivtop += row[i] * v[i]
            Pl[attr] = (self.bHash ** Mivtop) * (attrHash ** s_prime)
            Pr[attr] = mpk['g2'] ** s_prime
        CT1 = mpk['g2'] ** s
        CT2 = self._H2(mpk, t) ** s
        return {
            'policy': policy,
            'R': R,
            'Rt': Rt,
            'A': A,
            'Apr': Apr,
            'B': B,
            'C': C,
            'Ct': Ct,
            'Pl': Pl,
            'Pr': Pr,
            'CT1': CT1,
            'CT2': CT2,
            't': t
        }

    def tr_encrypt(self, mpk, M, policy_str, target_pos, t):
        i_target, j_target = target_pos
        policy = self.msp.createPolicy(policy_str)
        mono = self.msp.convert_policy_to_msp(policy)
        num_cols = self.msp.len_longest_row
        tau = self.group.random(ZR)
        eta = self.group.random(ZR)
        a = self.group.random(ZR)
        b = self.group.random(ZR)
        c = self.group.random(ZR)
        s_rows = [self.group.random(ZR) for _ in range(self.d)]
        w = {y: [self.group.random(ZR), self.group.random(ZR), self.group.random(ZR)] for y in range(1, self.d + 1)}
        v = [self.group.random(ZR) for _ in range(num_cols)]
        T=self.group.random(ZR)
        s = self.group.random(ZR)
        v[0] = s
        v1 = [a, self.group.init(ZR, 0), c]
        v2 = [self.group.init(ZR, 0), b, c]
        v3 = [-(b * c), -(a * c), (a * b)]
        v_c = [self.group.random(ZR), self.group.random(ZR), self.group.random(ZR)]
        R = {}
        Rt = {}
        A = {}
        Apr = {}
        B = {}
        for x in range(1, self.d + 1):
            qx = self.group.random(ZR)
            qx1 = self.group.random(ZR)
            qx2 = self.group.random(ZR)
            
            if x < i_target:
                vx=[self.group.random(ZR) for _ in range(3)]
                B[x]= self.group.pair_prod(mpk['g1'], mpk['g2']) ** T
            elif x == i_target:
                vx = [qx * v1[k]+qx2*v3[k] + qx1 * v2[k] for k in range(3)]
            elif x > i_target:
                vx = [qx * v1[k] + qx1 * v2[k] for k in range(3)]
               
            dot_v = vx[0] * v_c[0] + vx[1] * v_c[1] + vx[2] * v_c[2]
            if x>=i_target:
                B[x] = M * (mpk['G'][x - 1] ** (s_rows[x - 1] * tau * dot_v)) * (mpk['e_gg_alpha'] ** s)
            R[x] = [mpk['E'][x - 1] ** (s_rows[x - 1] * vx[k]) for k in range(3)]
            Rt[x] = [R[x][k] ** eta for k in range(3)]
            A[x] = mpk['g2'] ** (s_rows[x - 1] * tau * dot_v)
            Apr[x] =(mpk['g1'] ** (s_rows[x - 1] * tau * dot_v))

            
        C = {}
        Ct = {}
        
        for y in range(1, self.d + 1):
            vc4 = self.group.random(ZR)
            vc_prime = [v_c[k] + vc4 * v3[k] for k in range(3)]
            current_vc = vc_prime if y < j_target else v_c
            C[y] = [mpk['F'][y - 1] ** (tau * current_vc[k]) * (mpk['g1'] ** (w[y][k] * eta)) for k in range(3)]
            Ct[y] = [mpk['g1'] ** w[y][k] for k in range(3)]
        Pl = {}
        Pr = {}
        s_prime = self.group.random(ZR)
        for attr, row in mono.items():
            attr_stripped = self.msp.strip_index(attr)
            attrHash = self.group.hash(attr_stripped, G1)
            Mivtop = self.group.init(ZR, 0)
            for i in range(len(row)):
                Mivtop += row[i] * v[i]
            Pl[attr] = (self.bHash ** Mivtop) * (attrHash ** s_prime)
            Pr[attr] = mpk['g2'] ** s_prime
        CT1 = mpk['g2'] ** s
        CT2 = self._H2(mpk, t) ** s
        return {
            'policy': policy,
            'R': R,
            'Rt': Rt,
            'A': A,
            'Apr': Apr,
            'B': B,
            'C': C,
            'Ct': Ct,
            'Pl': Pl,
            'Pr': Pr,
            'CT1': CT1,
            'CT2': CT2,
            't': t
        }


    def trace(self, mpk, policy_str, t, decoder_oracle, W=20, epsilon=0.1):
        N = self.N
        p_estimates = []
        for k in range(1, N + 2):
            success_count = 0
            for _ in range(W):
                M_test = self.group.random(GT)
                target_pos = self._index_to_pos(k)
                if target_pos is None:
                    target_pos = (self.d + 1, 1)
                ct_probe = self.tr_encrypt(mpk, M_test, policy_str, target_pos, t)
                dec = None
                try:
                    dec = decoder_oracle(ct_probe)
                except Exception:
                    dec = None
                if dec == M_test:
                    success_count += 1
            p_estimates.append(success_count / float(W))
        traitors = []
        threshold = epsilon / (4 * N)
        for k in range(1, N + 1):
            gap = p_estimates[k - 1] - p_estimates[k]
            if gap >= threshold:
                target_pos = self._index_to_pos(k)
                if target_pos:
                    traitors.append(target_pos)
        return traitors

    def key_update(self, mpk, msk, t):
        cover_nodes = self.tree.cover_algorithm(self.revoked_users)
        uk = {}
        tr=self.group.init(ZR,t)
        for node in cover_nodes:
            theta = self.group.random(ZR)
            fn_t = node.polynomial_coeff * tr + msk['alpha']
            U1 = (mpk['g1'] ** fn_t) * (self._H2(mpk, t) ** theta)
            U2 = mpk['g2'] ** theta 
            uk[node.node_id] = {'U1': U1, 'U2': U2, 't': t}
        return uk

    def dec_key_gen(self, sk_user, uk_t):
        i, j = sk_user['pos']
        uid = self._user_id(i, j)
        path = self.tree.get_path(uid)
        node_id = None
        for node in path:
            if node.node_id in uk_t:
                node_id = node.node_id
                break
        if node_id is None:
            return None
        return {
            'sk1': sk_user['sk1'],
            'sk2': sk_user['sk2'],
            'sk3': sk_user['sk3'][node_id],
            'sk4_n': sk_user['sk4'][node_id],
            'sk5_n': sk_user['sk5'][node_id],
            'U1_n': uk_t[node_id]['U1'],
            'U2_n': uk_t[node_id]['U2'],
            'S': sk_user['S'],
            'pos': sk_user['pos'],
            't': uk_t[node_id]['t']
        }

    def decrypt(self, ct, dk):
        attrs = list(dk['S'])
        nodes = self.msp.prune(ct['policy'], attrs)
        if not nodes:
            return None
        coeffs = self.msp.getCoefficients(ct['policy'])
        num1=self.group.pair_prod(dk['sk4_n'], ct['CT1'])
        prod_rev_num = self.group.init(G1, 1)
        prod_rev_den = self.group.init(G1, 1)
        valid_attr_for_pr=None
        for node in nodes:
            attr = node.getAttributeAndIndex()
            attr_stripped = self.msp.strip_index(attr)
            if attr in ct['Pl'] and attr_stripped in dk['sk3'] and attr in ct['Pr']:
                if valid_attr_for_pr is None:
                    valid_attr_for_pr = attr
                num_shared = dk['sk3'][attr_stripped] ** coeffs[attr]
                
                den_shared = ct['Pl'][attr] ** coeffs[attr]
                prod_rev_num *= num_shared
                prod_rev_den *= den_shared
        pr1=self.group.pair_prod(prod_rev_num, ct['Pr'][valid_attr_for_pr])
        pr2=self.group.pair_prod(prod_rev_den,dk['sk5_n'])
        D_rev = num1*pr1 / pr2
        i, j = dk['pos']
        Rx = ct['R'][i]
        Rtx = ct['Rt'][i]
        Cy = ct['C'][j]
        Cty = ct['Ct'][j]
        num_tr = self.group.pair_prod(Rx[0], Cy[0]) * self.group.pair_prod(Rx[1], Cy[1]) * self.group.pair_prod(Rx[2], Cy[2])
        den_tr = self.group.pair_prod(Rtx[0], Cty[0]) * self.group.pair_prod(Rtx[1], Cty[1]) * self.group.pair_prod(Rtx[2], Cty[2])
        D_tr = num_tr / den_tr
        num_rev = self.group.pair_prod(dk['U1_n'], ct['CT1'])
        den_rev = self.group.pair_prod(ct['CT2'], dk['U2_n'])
        E_rev = num_rev / den_rev
       
        t_zr =self.group.init(ZR, dk['t'])
        if t_zr==self.group.init(ZR, 1):
            return None
        K_rev = (D_rev ** (t_zr / (t_zr - 1))) * (E_rev ** (1 / (1 - t_zr)))
        term1 = self.group.pair_prod(ct['Apr'][i], dk['sk2'])
        term2 = self.group.pair_prod(dk['sk1'], ct['A'][i])
        D_key = ct['B'][i] * (term1 / term2)
        D_final = D_key * D_tr
        M_rec = D_final / K_rev
        return M_rec
