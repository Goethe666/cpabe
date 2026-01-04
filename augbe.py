from charm.toolbox.pairinggroup import PairingGroup, G1, G2, GT, ZR
from charm.toolbox.secretutil import SecretUtil


class AugBE:
    """
    Augmented Broadcast Encryption (AugBE) - 非对称双线性群实现

    说明：
    - 使用 Type-3 配对（如 MNT224），配对为 e: G1 × G2 → GT
    - 参与配对的元素必须在正确的群域：左侧在 G1，右侧在 G2
    - 保持向量操作简单一致（二维向量，分量在 ZR）
    - 修正了论文中在非对称场景的若干域错配，使代数逻辑（尤其是 u_k 项约消）成立
    """

    def __init__(self, m=2, group_type='MNT224'):
        self.group = PairingGroup(group_type)
        self.m = m
        # Public parameters
        self.pk = None
        # Secret keys for all users (x,y)
        self.sk = {}

    # ---------------------- Utilities ----------------------
    def _zr(self, v):
        if isinstance(v, int):
            return self.group.init(ZR, v)
        return v

    def _vec_rand(self):
        return [self.group.random(ZR), self.group.random(ZR)]

    def _dot(self, a, b):
        return a[0] * b[0] + a[1] * b[1]

    def _pair_vec_prod(self, vec_g1, vec_g2):
        """
        逐分量配对并相乘：∏ e(vec_g1[k], vec_g2[k])。
        要求 vec_g1 为 G1 元素向量，vec_g2 为 G2 元素向量。
        """
        res = self.group.init(GT, 1)
        # 支持长度为 2 的向量（论文非对称方案的向量维度要求）
        for k in range(2):
            res *= self.group.pair_prod(vec_g1[k], vec_g2[k])
        return res

    # ---------------------- Setup ----------------------
    def setup(self):
        """
        Setup(m)：生成公钥与每个用户的密钥。
        用户按 (x, y) 编号，x ∈ [1..m], y ∈ [1..m]，总数 N = m^2。
        """
        g1 = self.group.random(G1)
        g2 = self.group.random(G2)
        e_gg = self.group.pair_prod(g1, g2)

        # 随机指数与公开参数生成
        r = [self.group.random(ZR) for _ in range(self.m)]         # r_x
        c = [self.group.random(ZR) for _ in range(self.m)]         # c_y
        alpha = [self.group.random(ZR) for _ in range(self.m)]     # alpha_x
        u = [self.group.random(G2) for _ in range(self.m)]         # u_y in G2

        E = [g1 ** r_x for r_x in r]                               # E_x = g1^{r_x}
        G = [e_gg ** a_x for a_x in alpha]                         # G_x = e(g1,g2)^{alpha_x}
        H = [g2 ** c_y for c_y in c]                               # H_y = g2^{c_y}

        self.pk = {
            'g1': g1,
            'g2': g2,
            'E': E,
            'G': G,
            'H': H,
            'u': u,
            'm': self.m
        }

        # 每个用户的密钥
        # 存储功能性组件以及解密所需的标量
        for x in range(1, self.m + 1):
            for y in range(1, self.m + 1):
                sigma_xy = self.group.random(ZR)
                # K_main 不直接用于解密，只是展示结构：g2^{alpha_x} * g2^{r_x c_y} * u_y^{sigma}
                K_main = (self.pk['g2'] ** alpha[x - 1]) * (self.pk['g2'] ** (r[x - 1] * c[y - 1])) * (u[y - 1] ** sigma_xy)
                # 修正：K_sigma 必须在 G1 中，以便与 G2 中的 T_x 配对（e(G1, G2)）
                K_sigma = self.pk['g1'] ** sigma_xy
                U_sigma = u[y - 1] ** sigma_xy

                self.sk[(x, y)] = {
                    'x': x,
                    'y': y,
                    'alpha_x': alpha[x - 1],
                    'r_x': r[x - 1],
                    'c_y': c[y - 1],
                    'u_y': u[y - 1],
                    'sigma': sigma_xy,
                    'K_main': K_main,
                    'K_sigma': K_sigma,
                    'U_sigma': U_sigma
                }

        return self.pk, self.sk

    # ---------------------- Encrypt ----------------------
    def encrypt(self, M, S_rows, i, j):
        """
        为接收者集合 S_rows 加密消息 M，目标索引为 (i, j)。

        输入：
          - M: GT 群元素（待加密）
          - S_rows: 字典 { x: [y1, y2, ...] }，表示每一行的接收者列索引
          - i, j: 目标位置（参考论文）

        输出：包含行与列组件的密文 Ctx。
        """
        g1 = self.pk['g1']
        g2 = self.pk['g2']
        E=self.pk['E']
        H=self.pk['H']
        G = self.pk['G']  # list over rows
        u = self.pk['u']  # list over columns

        # 随机标量与二维向量（保持论文 5.2 的二维结构）
        t = self.group.random(ZR)
        eta = self.group.random(ZR)
        s = [self.group.random(ZR) for _ in range(self.m)]
        w = [self._vec_rand() for _ in range(self.m)]  # w_y 为二维向量

        v1 = self._vec_rand()
        # Choose v2 orthogonal to v1: pick v2 = (v1[1], -v1[0])
        v2 = [v1[1], -v1[0]]
        v_c = self._vec_rand()
        v_cr = self.group.random(ZR)
        # v_c' = v_c + v_cr * v2
        v_c_prime = [v_c[0] + v_cr * v2[0], v_c[1] + v_cr * v2[1]]

        # 行组件（与行 x 相关）
        R = {}
        Rt = {}
        A = {}
        B = {}
        T = {}

        for x in range(1, self.m + 1):
            if x < i:
                # Random components
                z_x = self._vec_rand()
                a_x = self.group.random(ZR)
                b_x = self.group.random(ZR)
                # 向量指数：R_x, \tilde{R}_x 为 G1^2
                R[x] = [g1 ** z_x[0], g1 ** z_x[1]]
                Rt[x] = [R[x][0] ** eta, R[x][1] ** eta]
                A[x] = g1 ** a_x  # place in G1 to allow pairing with G2 on decrypt side
                # 修正：x<i 时 B_x 不包含消息 M
                B[x] = (G[x - 1] ** b_x)
                # T_x uses product over recipients in this row (S_x)
                Sx = S_rows.get(x, [])
                prod_u = self.group.init(G2, 1)
                for yk in Sx:
                    prod_u = prod_u * u[yk - 1]
                T[x] = prod_u ** a_x

            elif x == i:
                # Target row
                v_i = self._vec_rand()
                val = s[x - 1] * t * self._dot(v_i, v_c)
                # 向量指数：R_i = g1^{r_i s_i v_i}，通过 E_x = g1^{r_x} 的逐分量幂实现
                R[x] = [E[x - 1] ** (s[x - 1] * v_i[0]),
                        E[x - 1] ** (s[x - 1] * v_i[1])]
                Rt[x] = [R[x][0] ** eta, R[x][1] ** eta]
                A[x] = g1 ** val
                B[x] = M * (G[x - 1] ** val)
                Sx = S_rows.get(x, [])
                prod_u = self.group.init(G2, 1)
                for yk in Sx:
                    prod_u = prod_u * u[yk - 1]
                T[x] = prod_u ** val

            else:
                # Rows after target
                v_x_prime = self.group.random(ZR)  # 标量 v_x'（论文 5.2）
                v_x = [v1[0] * v_x_prime, v1[1] * v_x_prime]
                val = s[x - 1] * t * self._dot(v_x, v_c)
                R[x] = [E[x - 1] ** (s[x - 1] * v_x[0]),
                        E[x - 1] ** (s[x - 1] * v_x[1])]
                Rt[x] = [R[x][0] ** eta, R[x][1] ** eta]
                A[x] = g1 ** val
                B[x] = M * (G[x - 1] ** val)
                Sx = S_rows.get(x, [])
                prod_u = self.group.init(G2, 1)
                for yk in Sx:
                    prod_u = prod_u * u[yk - 1]
                T[x] = prod_u ** val

        # Column components
        C = {}
        Ct = {}
        for y in range(1, self.m + 1):
            vy = w[y - 1]  # 二维向量
            if y < j:
                base = v_c_prime
                C[y] = [ (H[y - 1] ** (t * base[0])) * (g2 ** (eta * vy[0])),
                         (H[y - 1] ** (t * base[1])) * (g2 ** (eta * vy[1])) ]
            else:
                base = v_c
                C[y] = [ (H[y - 1] ** (t * base[0])) * (g2 ** (eta * vy[0])),
                         (H[y - 1] ** (t * base[1])) * (g2 ** (eta * vy[1])) ]
            Ct[y] = [g2 ** vy[0], g2 ** vy[1]]

        Ctx = {
            'R': R,
            'Rt': Rt,
            'A': A,
            'B': B,
            'T': T,
            'C': C,
            'Ct': Ct,
            'S_rows': S_rows,
            'i': i,
            'j': j,
        }
        return Ctx

    # ---------------------- Decrypt ----------------------
    def decrypt(self, Ctx, user_xy):
        """
        为用户 (x, y) 解密密文 Ctx。
        返回恢复出的 GT 元素。
        """
        x, y = user_xy
        if (x, y) not in self.sk:
            return None
        sk = self.sk[(x, y)]
        g2 = self.pk['g2']

        # Collect row and column components
        R_x = Ctx['R'][x]
        Rt_x = Ctx['Rt'][x]
        A_x = Ctx['A'][x]    # in G1
        B_x = Ctx['B'][x]    # in GT
        T_x = Ctx['T'][x]    # in G2
        C_y = Ctx['C'][y]    # in G2
        Ct_y = Ctx['Ct'][y]  # in G2

        # 构造 K'_(x,y) = g2^{alpha_x} * g2^{r_x c_y} * Π_k (u_k^{sigma_{x,y}})
        Kprime_base = (self.pk['g2'] ** sk['alpha_x']) * (self.pk['g2'] ** (sk['r_x'] * sk['c_y']))
        prod_u_sigma = self.group.init(G2, 1)
        for yk in Ctx['S_rows'].get(x, []):
            prod_u_sigma = prod_u_sigma * (self.pk['u'][yk - 1] ** sk['sigma'])
        Kprime = Kprime_base * prod_u_sigma

        # 计算 M 的核心因子（保持配对域正确）：
        # factor1 = e(A_x, Kprime) / e(g1^{sigma_{x,y}}, T_x)
        num1 = self.group.pair_prod(A_x, Kprime)
        # 修正配对次序：K_sigma 在 G1，T_x 在 G2（从而 u_k 项可约消）
        den1 = self.group.pair_prod(sk['K_sigma'], T_x)
        factor1 = num1 / den1

        # 列相关因子：factor2 = ∏ e(R_x[k], C_y[k]) / ∏ e(Rt_x[k], Ct_y[k])
        num2 = self._pair_vec_prod(R_x, C_y)
        den2 = self._pair_vec_prod(Rt_x, Ct_y)
        factor2 = num2 / den2

        M_rec = B_x / factor1 * factor2
        return M_rec


def demo():
    group = PairingGroup('MNT224')
    scheme = AugBE(m=2, group_type='MNT224')
    pk, sk = scheme.setup()

    # Message to encrypt
    M = group.random(GT)

    # Recipients per row: deliver to row 1 users y>=1 and row 2 users y>=1
    S_rows = {1: [1, 2], 2: [1, 2]}
    i, j = 1, 1  # target index

    Ctx = scheme.encrypt(M, S_rows, i, j)

    # ---- 打印密文组件 ----
    print("=== 密文组件概览 ===")
    print(f"m={scheme.m}, i={i}, j={j}")
    print(f"S_rows={S_rows}")

    def _print_vec_g1(label, vec):
        if isinstance(vec, list):
            print(f"  {label} (len={len(vec)}):")
            for k, e in enumerate(vec):
                print(f"    [{k}] G1: {e}")
        else:
            print(f"  {label} (scalar G1): {vec}")

    def _print_vec_g2(label, vec):
        if isinstance(vec, list):
            print(f"  {label} (len={len(vec)}):")
            for k, e in enumerate(vec):
                print(f"    [{k}] G2: {e}")
        else:
            print(f"  {label} (scalar G2): {vec}")

    print("-- 行组件 --")
    for x in range(1, scheme.m + 1):
        print(f"Row x={x}:")
        _print_vec_g1("R[x]", Ctx['R'][x])
        _print_vec_g1("Rt[x]", Ctx['Rt'][x])
        print(f"  A[x] G1: {Ctx['A'][x]}")
        print(f"  B[x] GT: {Ctx['B'][x]}")
        _print_vec_g2("T[x]", Ctx['T'][x])

    print("-- 列组件 --")
    for y in range(1, scheme.m + 1):
        print(f"Col y={y}:")
        _print_vec_g2("C[y]", Ctx['C'][y])
        _print_vec_g2("Ct[y]", Ctx['Ct'][y])

    # Try to decrypt for user (x=1, y=1)
    M1 = scheme.decrypt(Ctx, (1, 1))
    print("Decrypt (1,1) success:", M1 == M)
    print("Recovered M1:", M1)
    print("Original  M :", M)

    # Try to decrypt for user (x=2, y=2)
    M2 = scheme.decrypt(Ctx, (2, 2))
    print("Decrypt (2,2) success:", M2 == M)
    print("Recovered M2:", M2)
def index_to_position(idx, m):
    """将一维索引 idx ∈ [1..m^2] 转为二维位置 (i', j')。"""
    i_prime = ((idx - 1) // m) + 1
    j_prime = ((idx - 1) % m) + 1
    return i_prime, j_prime

def augbe_success_rate_test(m=2, W=20):
    """
    AugBE 成功率测试：
    - S 设为全体用户集合；
    - 对 i = 1..N+1（N=m^2）以及 W 次试验，统计 D(ct)=M 的比例。
    说明：D 在此实现为合法用户 (1,1) 的解密器，以验证加密/解密正确性。
    """
    group = PairingGroup('MNT224')
    scheme = AugBE(m=m, group_type='MNT224')
    pk, sk = scheme.setup()

    # 全体用户集合 S_rows：每行包含所有列
    S_rows = {x: list(range(1, m + 1)) for x in range(1, m + 1)}
    N = m * m
    print("=== AugBE 成功率测试 ===")
    print(f"m={m}, N={N}, W={W}, S_rows={S_rows}")

    for idx in range(1, N + 2):  # 包含 N+1 次
        if idx <= N:
            i_prime, j_prime = index_to_position(idx, m)
        else:
            # 消息隐藏测试：映射到不存在的索引 (m+1, 1)
            i_prime, j_prime = m + 1, 1

        cnt = 0
        for _ in range(W):
            M = group.random(GT)
            ct = scheme.encrypt(M, S_rows, i_prime, j_prime)
            # 海盗解码器 D：这里用合法用户 (1,1) 的解密器代表
            M_rec = scheme.decrypt(ct, (1, 1))
            if M_rec == M:
                cnt += 1
        hat_p = cnt / float(W)
        print(f"i={idx} -> pos=({i_prime},{j_prime}), \u005Chat p_i = {hat_p:.3f}")


if __name__ == '__main__':
    demo()
 
    augbe_success_rate_test(m=2, W=20)

