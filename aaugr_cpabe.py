import random
from charm.toolbox.pairinggroup import PairingGroup, ZR, G1,G2, GT, pair
from MSP import MSP
class AugRCPABE:
    
    def __init__(self,m=2, group_name='MNT224'):
        self.group = PairingGroup(group_name)
        self.m = m
        self.N = self.m * self.m
        self.msp = MSP(self.group, verbose=False)

    def setup_A(self):
        m=self.m
        g1 = self.group.random(G1)
        g2 = self.group.random(G2)
        f = self.group.random(G1)
        h = self.group.random(G1)
        G = self.group.random(G1)
        H = self.group.random(G1)
        
        e_gg = pair(g1, g2)
        fx=[self.group.random(G1) for _ in range(m)]
        alpha = [self.group.random(ZR) for _ in range(m)]
        r = [self.group.random(ZR) for _ in range(m)]
        z = [self.group.random(ZR) for _ in range(m)]
        c = [self.group.random(ZR) for _ in range(m)]
        Ei = [e_gg ** alpha[i] for i in range(m)]
        Gi = [g1 ** r[i] for i in range(m)]
        Zi = [g1 ** z[i] for i in range(m)]
        Hj = [g2 ** c[j] for j in range(m)]
        pp = {
            'g1': g1,
            'g2': g2,   
            'f': f,
            'h': h,
            'fx': fx,
            'Ei': Ei,
            'Gi': Gi,
            'Zi': Zi,
            'Hj': Hj,
            'G': G,
            'H': H,
        }
        msk = {
            'alpha': alpha,
            'r': r,
            'c': c,
            'ctr': 0,
        }
        return pp, msk

    def _attr_to_zr(self, x):
        return self.group.hash(str(x), ZR)

    def _random_vec(self):
        return [self.group.random(ZR) for _ in range(3)]

    def _dot(self, v1, v2):
        return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]

    def _ctr_to_pos(self, ctr):
        i = (ctr // self.m) + 1
        j = (ctr % self.m) + 1
        return (i, j)

    def _index_to_pos(self, k):
        if 1 <= k <= self.N:
            return self._ctr_to_pos(k - 1)
        return None

    def KeyGen_A(self, pp, msk, S):
        msk['ctr'] += 1
        idx = (msk['ctr'] - 1) % (self.m * self.m)
        i, j = self._ctr_to_pos(idx)
        sigma = self.group.random(ZR)
        K_ij = (pp['g1'] ** msk['alpha'][i - 1]) * ((pp['Gi'][i - 1] ** msk['c'][j - 1])) * ((pp['f'] * pp['fx'][j - 1]) ** sigma)
        Kp_ij = pp['g2'] ** sigma
        Kpp_ij = pp['Zi'][i - 1] ** sigma
        barK = {}
        for jp in range(1, self.m + 1):
            if jp != j:
                barK[jp] = pp['fx'][jp - 1] ** sigma
        K_attr = {}
        Kp_attr = {}
        for x in S:
            delta = self.group.random(ZR)
            K_attr[x] = pp['g2'] ** delta
            Hz = (pp['H'] ** self._attr_to_zr(x)) * pp['h']
            Kp_attr[x] = (Hz ** delta) * (pp['G'] ** (-sigma))
        return {
            'i': i,
            'j': j,
            'sigma': sigma,
            'S': set(S),
            'K_ij': K_ij,
            'Kp_ij': Kp_ij,
            'Kpp_ij': Kpp_ij,
            'barK': barK,
            'K_attr': K_attr,
            'Kp_attr': Kp_attr,
        }

    def Encrypt_A(self, pp, M, R, policy_str, target_pos=(1,1)):
        if target_pos is None:
            target_pos=(self.m+1,1)
        policy = self.msp.createPolicy(policy_str)
        mono = self.msp.convert_policy_to_msp(policy)
        n = self.msp.len_longest_row
        pi = self.group.random(ZR)
        u = [pi] + [self.group.random(ZR) for _ in range(n - 1)]
        xi = {attr: self.group.random(ZR) for attr in mono.keys()}
        tau = self.group.random(ZR)
        kappa = self.group.random(ZR)
        s=[self.group.random(ZR) for _ in range(1,self.m+1)]
        t=[self.group.random(ZR) for _ in range(1,self.m+1)]
        a = self.group.random(ZR)
        b = self.group.random(ZR)
        c = self.group.random(ZR)
        v1 = [a, self.group.init(ZR, 0), c]
        v2 = [self.group.init(ZR, 0), b, c]
        v3 = [v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0]]
        vc = self._random_vec()
        w = {j: self._random_vec() for j in range(1, self.m + 1)}
        P = {}
        Pp = {}
        Ppp = {}
        for attr, row in mono.items():
            acc = self.group.init(ZR, 0)
            for j in range(len(row)):
                acc += self.group.init(ZR, int(row[j])) * u[j]
            P[attr] = (pp['f'] ** acc) * (pp['G'] ** xi[attr])
            rho_z = self._attr_to_zr(self.msp.strip_index(attr))
            Pp[attr] = ((pp['H'] ** rho_z) * pp['h']) ** (-xi[attr])
            Ppp[attr] = pp['g2'] ** xi[attr]
        revoked_set = set(R or [])
        Ri = {}
        Rp = {}
        Q={}
        Qp={}
        Qpp={}
        T={}
        
        for i in range(1, self.m + 1):
            bar_R_i = [j for j in range(1, self.m + 1) if (i, j) not in revoked_set]
            F_prod = pp['f']
            for jp in bar_R_i:
                F_prod = F_prod * pp['fx'][jp - 1]
            qx = self.group.random(ZR)
            qx1 = self.group.random(ZR)
            if i<target_pos[0]:
                vi=self._random_vec()
                s_2=self.group.random(ZR)
                Ri[i] = [pp['g1'] ** vi[k] for k in range(3)]
                Rp[i] = [Ri[i][k] ** kappa for k in range(3)]
                Q[i] = pp['g2'] ** (s[i - 1])
                Qp[i] = (F_prod)**(s[i - 1])*pp['Zi'][i - 1]**t[i - 1]*pp['f']**pi
                Qpp[i] = pp['g2'] ** (t[i - 1])
                T[i] = (pp['Ei'][i - 1] ** s_2)
            elif i==target_pos[0]:
                vi=self._random_vec()
                dot_vi_vc = self._dot(vi, vc)
                Ri[i] = [pp['Gi'][i - 1] ** (s[i - 1] * vi[k]) for k in range(3)]    
                Rp[i] = [Ri[i][k] ** kappa for k in range(3)]
                Q[i] = pp['g2'] ** (tau * s[i - 1] * dot_vi_vc)
                Qp[i] = (F_prod)**(tau*s[i - 1]*dot_vi_vc)*pp['Zi'][i - 1]**t[i - 1]*pp['f']**pi
                Qpp[i] = pp['g2'] ** (t[i - 1])
                T[i] = M * (pp['Ei'][i - 1] ** (tau * s[i - 1] * dot_vi_vc))
            elif i>target_pos[0]:
                vi=[qx*v1[k]+qx1*v2[k] for k in range(3)]
                dot_vi_vc = self._dot(vi, vc)
                Ri[i] = [pp['Gi'][i - 1] ** (s[i - 1] * vi[k]) for k in range(3)]    
                Rp[i] = [Ri[i][k] ** kappa for k in range(3)]
                Q[i] = pp['g2'] ** (tau * s[i - 1] * dot_vi_vc)
                Qp[i] = (F_prod)**(tau*s[i - 1]*dot_vi_vc)*pp['Zi'][i - 1]**t[i - 1]*pp['f']**pi
                Qpp[i] = pp['g2'] ** (t[i - 1])
                T[i] = M * (pp['Ei'][i - 1] ** (tau * s[i - 1] * dot_vi_vc))

        C={}
        Cp={}
        for j in range(1, self.m + 1):
            miu=self.group.random(ZR)
            vc_prime = [vc[k] + miu * v3[k] for k in range(3)] 
            current_vc = vc_prime if j < target_pos[1] else vc
            C[j] = [pp['Hj'][j - 1] ** (tau * current_vc[k]) * (pp['g2'] ** (w[j][k] * kappa)) for k in range(3)]
            Cp[j] = [pp['g2'] ** w[j][k] for k in range(3)]
        return {
            'policy': policy,
            'Ri':Ri,
            'Rp':Rp,
            'Q':Q,
            'Qp':Qp,
            'Qpp':Qpp,
            'T':T,
            'C':C,
            'Cp':Cp,
            'P':P,
            'Pp':Pp,
            'Ppp':Ppp,
            'revoked': revoked_set,
            'pos': target_pos,
        }

    def Decrypt_A(self, pp, CT, SK):
        i=SK['i']
        j=SK['j']
        attrs = list(SK['S'])
        pruned = self.msp.prune(CT['policy'], attrs)
        if not pruned:
            return None
        coeffs = self.msp.getCoefficients(CT['policy'])
        Dp_num = self.group.init(GT, 1)
        for node in pruned:
            aidx = node.getAttributeAndIndex()
            aplain = self.msp.strip_index(aidx)
            w = coeffs[aidx]
            Dp_num *= (pair(SK['Kp_ij'], CT['P'][aidx]) * pair(SK['K_attr'][aplain], CT['Pp'][aidx]) * pair(SK['Kp_attr'][aplain], CT['Ppp'][aidx])) ** w
        bar_R_i = set([jp for jp in range(1, self.m + 1) if (i, jp) not in CT.get('revoked', set())])
        barK_i = SK['K_ij']
        for jp in bar_R_i:
            if jp != j:
                barK_i = barK_i * SK['barK'].get(jp, pp['g1'] ** self.group.init(ZR, 0))
        Di_pre = (pair(barK_i, CT['Q'][i]) * pair(SK['Kpp_ij'], CT['Qpp'][i])) / pair(SK['Kp_ij'], CT['Qp'][i])
        num_pair = self.group.init(GT, 1)
        den_pair = self.group.init(GT, 1)
        for k in range(3):
            num_pair *= pair(CT['Rp'][i][k], CT['Cp'][j][k])
            den_pair *= pair(CT['Ri'][i][k], CT['C'][j][k])
        Di = Di_pre * (num_pair / den_pair)
        denom = Dp_num * Di
        Ti = CT['T'][i]
        return Ti / denom

    def Trace(self, pp, R_D, policy_str, epsilon, decoder_oracle, W=20):
        N = self.N
        p_estimates = []
        for k in range(1, N + 2):
            success_count = 0
            target_pos = self._index_to_pos(k)
            for _ in range(W):
                M_test = self.group.random(GT)
                ct_probe = self.Encrypt_A(pp, M_test, R_D, policy_str, target_pos)
                dec = decoder_oracle(ct_probe)
                if dec == M_test:
                    success_count += 1
            p_estimates.append(success_count / float(W))
        traitors = []
        threshold = epsilon / (4 * N) if N > 0 else epsilon
        for k in range(1, N + 1):
            gap = p_estimates[k - 1] - p_estimates[k]
            if gap >= threshold:
                traitors.append(k)
        return traitors
