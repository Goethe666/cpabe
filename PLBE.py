from charm.toolbox.pairinggroup import PairingGroup, G1, G2, GT, ZR


class PLBE:
    def __init__(self, m=2, group_type='SS512'):
        self.group = PairingGroup(group_type)
        self.m = m
        self.pk = None
        self.sk = {}

    def _vec3(self):
        return [self.group.random(ZR), self.group.random(ZR), self.group.random(ZR)]

    def _dot3(self, a, b):
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

    def _pair_vec3(self, vec_g1, vec_g2):
        res = self.group.init(GT, 1)
        for k in range(3):
            res *= self.group.pair_prod(vec_g1[k], vec_g2[k])
        return res

    def setup_plbe(self):
        g1 = self.group.random(G1)
        g2 = self.group.random(G2)
        e_g1g2 = self.group.pair_prod(g1, g2)
        r = [self.group.random(ZR) for _ in range(self.m)]
        c = [self.group.random(ZR) for _ in range(self.m)]
        alpha = [self.group.random(ZR) for _ in range(self.m)]
        E = [g1 ** r_x for r_x in r]
        G = [e_g1g2 ** a_x for a_x in alpha]
        H = [g2 ** c_y for c_y in c]
        self.pk = {'g1': g1, 'g2': g2, 'E': E, 'G': G, 'H': H}
        for x in range(1, self.m + 1):
            for y in range(1, self.m + 1):
                Kxy = (g2 ** alpha[x - 1]) * (g2 ** (r[x - 1] * c[y - 1]))
                self.sk[(x, y)] = {'K': Kxy}
        return self.pk, self.sk

    def encrypt_plbe(self, M):
        if self.pk is None:
            self.setup_plbe()
        t = self.group.random(ZR)
        eta = self.group.random(ZR)
        s = [self.group.random(ZR) for _ in range(self.m)]
        w = {y: self._vec3() for y in range(1, self.m + 1)}
        a = self.group.random(ZR)
        b = self.group.random(ZR)
        c = self.group.random(ZR)
        v1 = [a, self.group.init(ZR, 0), c]
        v2 = [self.group.init(ZR, 0), b, c]
        v3 = [-(b * c), -(a * c), (a * b)]
        vc = self._vec3()
        R = {}
        Rt = {}
        A = {}
        B = {}
        for x in range(1, self.m + 1):
            q1 = self.group.random(ZR)
            q2 = self.group.random(ZR)
            vx = [q1 * v1[0] + q2 * v2[0], q1 * v1[1] + q2 * v2[1], q1 * v1[2] + q2 * v2[2]]
            dot_v = self._dot3(vx, vc)
            R[x] = [self.pk['E'][x - 1] ** (s[x - 1] * vx[0]), self.pk['E'][x - 1] ** (s[x - 1] * vx[1]), self.pk['E'][x - 1] ** (s[x - 1] * vx[2])]
            Rt[x] = [R[x][0] ** eta, R[x][1] ** eta, R[x][2] ** eta]
            A[x] = self.pk['g1'] ** (s[x - 1] * t * dot_v)
            B[x] = M * (self.pk['G'][x - 1] ** (s[x - 1] * t * dot_v))
        C = {}
        Ct = {}
        for y in range(1, self.m + 1):
            C[y] = [self.pk['H'][y - 1] ** (t * vc[0]) * (self.pk['g2'] ** (eta * w[y][0])), self.pk['H'][y - 1] ** (t * vc[1]) * (self.pk['g2'] ** (eta * w[y][1])), self.pk['H'][y - 1] ** (t * vc[2]) * (self.pk['g2'] ** (eta * w[y][2]))]
            Ct[y] = [self.pk['g2'] ** w[y][0], self.pk['g2'] ** w[y][1], self.pk['g2'] ** w[y][2]]
        return {'R': R, 'Rt': Rt, 'A': A, 'B': B, 'C': C, 'Ct': Ct}

    def decrypt_plbe(self, ct, x, y):
        if self.pk is None:
            self.setup_plbe()
        K = self.sk[(x, y)]['K']
        A_x = ct['A'][x]
        B_x = ct['B'][x]
        R_x = ct['R'][x]
        Rt_x = ct['Rt'][x]
        C_y = ct['C'][y]
        Ct_y = ct['Ct'][y]
        num1 = self.group.pair_prod(A_x, K)
        num2 = self._pair_vec3(R_x, C_y)
        den2 = self._pair_vec3(Rt_x,Ct_y)
        return B_x / num1 * (num2 / den2)

    def tr_encrypt_plbe(self, M, i, j):
        if self.pk is None:
            self.setup_plbe()
        t = self.group.random(ZR)
        eta = self.group.random(ZR)
        s = [self.group.random(ZR) for _ in range(self.m)]
        w = {y: self._vec3() for y in range(1, self.m + 1)}
        a = self.group.random(ZR)
        b = self.group.random(ZR)
        c = self.group.random(ZR)
        v1 = [a, self.group.init(ZR, 0), c]
        v2 = [self.group.init(ZR, 0), b, c]
        v3 = [-(b * c), -(a * c), (a * b)]
        vc = self._vec3()
        vc4 = self.group.random(ZR)
        vc_prime = [vc[0] + vc4 * v3[0], vc[1] + vc4 * v3[1], vc[2] + vc4 * v3[2]]
        R = {}
        Rt = {}
        A = {}
        B = {}
        for x in range(1, self.m + 1):
            if x < i:
                z1 = self._vec3()
                z2 = self.group.random(ZR)
                z3 = self.group.random(ZR)
                R[x] = [self.pk['g1'] ** z1[0], self.pk['g1'] ** z1[1], self.pk['g1'] ** z1[2]]
                Rt[x] = [R[x][0] ** eta, R[x][1] ** eta, R[x][2] ** eta]
                A[x] = self.pk['g1'] ** z2
                e_gg = self.group.pair_prod(self.pk['g1'], self.pk['g2'])
                B[x] = e_gg ** z3
            elif x == i:
                q1 = self.group.random(ZR)
                q2 = self.group.random(ZR)
                q3 = self.group.random(ZR)
                vx = [q1 * v1[0] + q2 * v2[0] + q3 * v3[0], q1 * v1[1] + q2 * v2[1] + q3 * v3[1], q1 * v1[2] + q2 * v2[2] + q3 * v3[2]]
                dot_v = self._dot3(vx, vc)
                R[x] = [self.pk['E'][x - 1] ** (s[x - 1] * vx[0]), self.pk['E'][x - 1] ** (s[x - 1] * vx[1]), self.pk['E'][x - 1] ** (s[x - 1] * vx[2])]
                Rt[x] = [R[x][0] ** eta, R[x][1] ** eta, R[x][2] ** eta]
                A[x] = self.pk['g1'] ** (s[x - 1] * t * dot_v)
                B[x] = M * (self.pk['G'][x - 1] ** (s[x - 1] * t * dot_v))
            else:
                q1 = self.group.random(ZR)
                q2 = self.group.random(ZR)
                vx = [q1 * v1[0] + q2 * v2[0], q1 * v1[1] + q2 * v2[1], q1 * v1[2] + q2 * v2[2]]
                dot_v = self._dot3(vx, vc)
                R[x] = [self.pk['E'][x - 1] ** (s[x - 1] * vx[0]), self.pk['E'][x - 1] ** (s[x - 1] * vx[1]), self.pk['E'][x - 1] ** (s[x - 1] * vx[2])]
                Rt[x] = [R[x][0] ** eta, R[x][1] ** eta, R[x][2] ** eta]
                A[x] = self.pk['g1'] ** (s[x - 1] * t * dot_v)
                B[x] = M * (self.pk['G'][x - 1] ** (s[x - 1] * t * dot_v))
        C = {}
        Ct = {}
        for y in range(1, self.m + 1):
            base = vc_prime if y < j else vc
            C[y] = [self.pk['H'][y - 1] ** (t * base[0]) * (self.pk['g2'] ** (eta * w[y][0])), self.pk['H'][y - 1] ** (t * base[1]) * (self.pk['g2'] ** (eta * w[y][1])), self.pk['H'][y - 1] ** (t * base[2]) * (self.pk['g2'] ** (eta * w[y][2]))]
            Ct[y] = [self.pk['g2'] ** w[y][0], self.pk['g2'] ** w[y][1], self.pk['g2'] ** w[y][2]]
        return {'R': R, 'Rt': Rt, 'A': A, 'B': B, 'C': C, 'Ct': Ct, 'i': i, 'j': j}

    def _index_to_pos(self, k):
        if k < 1 or k > self.m * self.m:
            return None
        i_prime = ((k - 1) // self.m) + 1
        j_prime = ((k - 1) % self.m) + 1
        return (i_prime, j_prime)

    def trace(self, decoder_oracle, W=20, epsilon=0.1):
        N = self.m * self.m
        p_estimates = []
        for k in range(1, N + 2):
            success_count = 0
            for _ in range(W):
                M_test = self.group.random(GT)
                pos = self._index_to_pos(k)
                if pos is None:
                    pos = (self.m + 1, 1)
                ct_probe = self.tr_encrypt_plbe(M_test, pos[0], pos[1])
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
                pos = self._index_to_pos(k)
                if pos:
                    traitors.append(pos)
        return traitors


    
