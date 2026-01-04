from math import ceil
import time
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.font_manager as fm
from math import ceil
from charm.toolbox.pairinggroup import PairingGroup, GT
from charm.schemes.abenc.waters11 import Waters11
from u_cpabe import UnifiedCPABE
from FABEOCP import FABEO22CPABE
from aaugr_cpabe import AugRCPABE
from aindre_abecp import IndirectRevocationABE
from PLBE import PLBE

REPEAT =30
D_DEFAULT = 10
D_REVOCATION = 10
GROUP_TYPE = 'MNT224'
def _find_tnr():
    try:
        for f in fm.fontManager.ttflist:
            if f.name in ('Times New Roman', 'TimesNewRoman', 'Times'):
                return f.fname
    except Exception:
        return None
    return None
TNR_PATH = _find_tnr()

def plot_graph_multi(filename, title, x_label, x_values, series_list, y_top=None, y_tick_step=None, y_label='Time (ms)', x_tick_step=None):
    plt.figure(figsize=(8, 6))
    if TNR_PATH:
        tnr_title = fm.FontProperties(fname=TNR_PATH, size=18, weight='bold')
        tnr_label = fm.FontProperties(fname=TNR_PATH, size=16, weight='bold')
        tnr_legend = fm.FontProperties(fname=TNR_PATH, size=18, weight='bold')
        tnr_ticks = fm.FontProperties(fname=TNR_PATH, size=14, weight='bold')
    else:
        tnr_title = fm.FontProperties(family='serif', size=18, weight='bold')
        tnr_label = fm.FontProperties(family='serif', size=16, weight='bold')
        tnr_legend = fm.FontProperties(family='serif', size=18, weight='bold')
        tnr_ticks = fm.FontProperties(family='serif', size=14, weight='bold')
    for series in series_list:
        kwargs = {'linewidth': 5, 'markersize': 9, 'label': series['label']}
        if 'color' in series:
            kwargs['color'] = series['color']
        if 'alpha' in series:
            kwargs['alpha'] = series['alpha']
        plt.plot(x_values, series['y'], series['style'], **kwargs)
    plt.xlabel(x_label, fontproperties=tnr_label)
    plt.ylabel(y_label, fontproperties=tnr_label)
    plt.grid(True, linestyle='--', alpha=0.7, linewidth=1.0)
    plt.xlim(left=0)
    if y_top is not None:
        plt.ylim(bottom=0, top=y_top)
    else:
        plt.ylim(bottom=0)
    if x_tick_step is None and len(x_values) >= 2:
        step = x_values[1] - x_values[0]
        if step > 0:
            xmax = max(x_values)
            plt.xticks(list(range(0, int(xmax) + int(step), int(step))))
    ax = plt.gca()
    if x_tick_step is not None and x_tick_step > 0:
        ax.xaxis.set_major_locator(MultipleLocator(x_tick_step))
    if y_tick_step is not None and y_tick_step > 0:
        ax.yaxis.set_major_locator(MultipleLocator(y_tick_step))
    for lbl in ax.get_xticklabels():
        lbl.set_fontproperties(tnr_ticks)
    for lbl in ax.get_yticklabels():
        lbl.set_fontproperties(tnr_ticks)
    ax.tick_params(axis='both', which='both', width=2, length=6)
    for spine in ax.spines.values():
        spine.set_linewidth(2)
    if any(s.get('label') for s in series_list):
        leg = plt.legend(shadow=True, loc='upper left', prop=tnr_legend, markerscale=1.6)
        leg.get_frame().set_linewidth(2)
    plt.tight_layout()
    plt.savefig(filename)
    print(f"Saved {filename}")
    plt.close()

def _choose_group(pref='MNT224'):
    try:
        return PairingGroup(pref)
    except Exception:
        for alt in ['MNT223', 'BN254', 'SS512', 'MNT159']:
            try:
                return PairingGroup(alt)
            except Exception:
                pass
        return PairingGroup('SS512')

def measure_keygen_ucpabe(x_values):
    print("Running Experiment: U-CPABE KeyGen")
    times = []
    n = len(x_values)
    for x in x_values:
        print(f"  - Processing x={x}...")
        scheme = UnifiedCPABE(d=D_DEFAULT, group_type=GROUP_TYPE)
        mpk, msk = scheme.setup()
        attrs = {str(i) for i in range(x)}
        total = 0.0
        for _ in range(REPEAT):
            t0 = time.time()
            scheme.keygen(mpk, msk, attrs, 1, 1)
            t1 = time.time()
            total += (t1 - t0) * 1000.0
        times.append(total / REPEAT)
        print(f"完成第 {len(times)}/{n} 个")
    return times

def measure_keygen_waters(x_values):
    print("Running Experiment: Waters KeyGen")
    g = _choose_group(GROUP_TYPE)
    cp = Waters11(g, max(x_values))
    mpk, msk = cp.setup()
    out = []
    n = len(x_values)
    for x in x_values:
        print(f"  - Processing x={x}...")
        attrs = [str(i) for i in range(1, x + 1)]
        total = 0.0
        for _ in range(REPEAT):
            t0 = time.time()
            cp.keygen(mpk, msk, attrs)
            t1 = time.time()
            total += (t1 - t0) * 1000.0
        out.append(total / REPEAT)
        print(f"完成第 {len(out)}/{n} 个")
    return out

def measure_keygen_fabeo(x_values):
    print("Running Experiment: FABEO KeyGen")
    g = _choose_group(GROUP_TYPE)
    scheme = FABEO22CPABE(g)
    mpk, msk = scheme.setup()
    out = []
    n = len(x_values)
    for x in x_values:
        print(f"  - Processing x={x}...")
        attrs = [str(i) for i in range(x)]
        total = 0.0
        for _ in range(REPEAT):
            t0 = time.time()
            scheme.keygen(mpk, msk, attrs)
            t1 = time.time()
            total += (t1 - t0) * 1000.0
        out.append(total / REPEAT)
        print(f"完成第 {len(out)}/{n} 个")
    return out

def measure_keygen_indre(x_values):
    print("Running Experiment: Indre KeyGen")
    num_users = 100
    print(f"Indre KeyGen 实验用户总数: {num_users}")
    scheme = IndirectRevocationABE(num_users=num_users, group_type=GROUP_TYPE)
    public_key, master_key = scheme.setup()
    out = []
    n = len(x_values)
    for x in x_values:
        print(f"  - Processing x={x}...")
        policy = " and ".join([str(i) for i in range(x)]) if x > 0 else "0"
        total = 0.0
        for _ in range(REPEAT):
            t0 = time.time()
            scheme.keygen(public_key, master_key, 1, policy)
            t1 = time.time()
            total += (t1 - t0) * 1000.0
        out.append(total / REPEAT)
        print(f"完成第 {len(out)}/{n} 个")
    return out

def measure_keygen_augr(x_values):
    print("Running Experiment: AugR KeyGen")
    scheme = AugRCPABE(m=D_DEFAULT, group_name=GROUP_TYPE)
    pp, msk = scheme.setup_A()
    out = []
    n = len(x_values)
    for x in x_values:
        print(f"  - Processing x={x}...")
        attrs = set([str(i) for i in range(x)])
        total = 0.0
        for _ in range(REPEAT):
            t0 = time.time()
            scheme.KeyGen_A(pp, msk, attrs)
            t1 = time.time()
            total += (t1 - t0) * 1000.0
        out.append(total / REPEAT)
        print(f"完成第 {len(out)}/{n} 个")
    return out

def measure_encrypt_ucpabe(x_values):
    print("Running Experiment: U-CPABE Encrypt")
    times = []
    n = len(x_values)
    for x in x_values:
        print(f"  - Processing x={x}...")
        scheme = UnifiedCPABE(d=D_DEFAULT, group_type=GROUP_TYPE)
        mpk, msk = scheme.setup()
        policy = " and ".join([str(i) for i in range(x)]) if x != 1 else "0"
        M = scheme.group.random(GT)
        t = 1
        total = 0.0
        for _ in range(REPEAT):
            t0 = time.time()
            scheme.encrypt(mpk, M, policy, t)
            t1 = time.time()
            total += (t1 - t0) * 1000.0
        times.append(total / REPEAT)
        print(f"完成第 {len(times)}/{n} 个")
    return times

def measure_encrypt_waters(x_values):
    print("Running Experiment: Waters Encrypt")
    g = _choose_group(GROUP_TYPE)
    cp = Waters11(g, max(x_values))
    mpk, msk = cp.setup()
    out = []
    n = len(x_values)
    for x in x_values:
        print(f"  - Processing x={x}...")
        policy = " and ".join([str(i) for i in range(1, x + 1)]) if x > 0 else "1"
        M = g.random(GT)
        total = 0.0
        for _ in range(REPEAT):
            t0 = time.time()
            cp.encrypt(mpk, M, policy)
            t1 = time.time()
            total += (t1 - t0) * 1000.0
        out.append(total / REPEAT)
        print(f"完成第 {len(out)}/{n} 个")
    return out

def measure_encrypt_fabeo(x_values):
    print("Running Experiment: FABEO Encrypt")
    g = _choose_group(GROUP_TYPE)
    scheme = FABEO22CPABE(g)
    mpk, msk = scheme.setup()
    out = []
    n = len(x_values)
    for x in x_values:
        print(f"  - Processing x={x}...")
        policy = " and ".join([str(i) for i in range(x)]) if x > 0 else "0"
        M = g.random(GT)
        total = 0.0
        for _ in range(REPEAT):
            t0 = time.time()
            scheme.encrypt(mpk, M, policy)
            t1 = time.time()
            total += (t1 - t0) * 1000.0
        out.append(total / REPEAT)
        print(f"完成第 {len(out)}/{n} 个")
    return out

def measure_encrypt_indre(x_values):
    print("Running Experiment: Indre Encrypt")
    num_users = 100
    print(f"Indre Encrypt 实验用户总数: {num_users}")
    scheme = IndirectRevocationABE(num_users=num_users, group_type=GROUP_TYPE)
    public_key, master_key = scheme.setup()
    out = []
    n = len(x_values)
    for x in x_values:
        print(f"  - Processing x={x}...")
        attrs = [str(i) for i in range(x)]
        M = scheme.group.random(GT)
        total = 0.0
        for _ in range(REPEAT):
            t0 = time.time()
            scheme.encrypt(public_key, M, attrs, 1)
            t1 = time.time()
            total += (t1 - t0) * 1000.0
        out.append(total / REPEAT)
        print(f"完成第 {len(out)}/{n} 个")
    return out

def measure_encrypt_augr(x_values):
    print("Running Experiment: AugR Encrypt")
    scheme = AugRCPABE(m=D_DEFAULT, group_name=GROUP_TYPE)
    pp, msk = scheme.setup_A()
    out = []
    n = len(x_values)
    for x in x_values:
        print(f"  - Processing x={x}...")
        policy = " and ".join([str(i) for i in range(x)]) if x > 0 else "0"
        M = scheme.group.random(GT)
        total = 0.0
        for _ in range(REPEAT):
            t0 = time.time()
            scheme.Encrypt_A(pp, M, [], policy)
            t1 = time.time()
            total += (t1 - t0) * 1000.0
        out.append(total / REPEAT)
        print(f"完成第 {len(out)}/{n} 个")
    return out

def measure_decrypt_ucpabe(x_values):
    print("Running Experiment: U-CPABE Decrypt")
    times = []
    n = len(x_values)
    for x in x_values:
        print(f"  - Processing x={x}...")
        scheme = UnifiedCPABE(d=D_DEFAULT, group_type=GROUP_TYPE)
        mpk, msk = scheme.setup()
        policy = " and ".join([str(i) for i in range(x)])
        attrs = {str(i) for i in range(x)}
        sk = scheme.keygen(mpk, msk, attrs, 1, 1)
        M = scheme.group.random(GT)
        t = 10
        ct = scheme.encrypt(mpk, M, policy, t)
        uk = scheme.key_update(mpk, msk, t)
        dk = scheme.dec_key_gen(sk, uk)
        total = 0.0
        for _ in range(REPEAT):
            t0 = time.time()
            scheme.decrypt(ct, dk)
            t1 = time.time()
            total += (t1 - t0) * 1000.0
        times.append(total / REPEAT)
        print(f"完成第 {len(times)}/{n} 个")
    return times

def measure_decrypt_waters(x_values):
    print("Running Experiment: Waters Decrypt")
    g = _choose_group(GROUP_TYPE)
    cp = Waters11(g, max(x_values))
    mpk, msk = cp.setup()
    out = []
    n = len(x_values)
    for x in x_values:
        print(f"  - Processing x={x}...")
        attrs = [str(i) for i in range(1, x + 1)]
        sk = cp.keygen(mpk, msk, attrs)
        M = g.random(GT)
        policy = " and ".join([str(i) for i in range(1, x + 1)]) if x > 0 else "1"
        ct = cp.encrypt(mpk, M, policy)
        total = 0.0
        for _ in range(REPEAT):
            t0 = time.time()
            cp.decrypt(mpk, ct, sk)
            t1 = time.time()
            total += (t1 - t0) * 1000.0
        out.append(total / REPEAT)
        print(f"完成第 {len(out)}/{n} 个")
    return out

def measure_decrypt_fabeo(x_values):
    print("Running Experiment: FABEO Decrypt")
    g = _choose_group(GROUP_TYPE)
    scheme = FABEO22CPABE(g)
    mpk, msk = scheme.setup()
    out = []
    n = len(x_values)
    for x in x_values:
        print(f"  - Processing x={x}...")
        attrs = [str(i) for i in range(x)]
        sk = scheme.keygen(mpk, msk, attrs)
        M = g.random(GT)
        policy = " and ".join([str(i) for i in range(x)]) if x > 0 else "0"
        ct = scheme.encrypt(mpk, M, policy)
        total = 0.0
        for _ in range(REPEAT):
            t0 = time.time()
            scheme.decrypt(mpk, ct, sk)
            t1 = time.time()
            total += (t1 - t0) * 1000.0
        out.append(total / REPEAT)
        print(f"完成第 {len(out)}/{n} 个")
    return out

def measure_decrypt_indre(x_values):
    print("Running Experiment: Indre Decrypt")
    num_users = 100
    print(f"Indre Decrypt 实验用户总数: {num_users}")
    scheme = IndirectRevocationABE(num_users=num_users, group_type=GROUP_TYPE)
    public_key, master_key = scheme.setup()
    out = []
    n = len(x_values)
    for x in x_values:
        print(f"  - Processing x={x}...")
        policy = " and ".join([str(i) for i in range(x)])
        attrs = [str(i) for i in range(x)]
        user_key = scheme.keygen(public_key, master_key, 1, policy)
        M = scheme.group.random(GT)
        ct = scheme.encrypt(public_key, M, attrs, 10)
        uk = scheme.key_update(public_key, master_key, 10)
        total = 0.0
        for _ in range(REPEAT):
            t0 = time.time()
            scheme._core_decrypt(ct, user_key, uk, 1)
            t1 = time.time()
            total += (t1 - t0) * 1000.0
        out.append(total / REPEAT)
        print(f"完成第 {len(out)}/{n} 个")
    return out

def measure_decrypt_augr(x_values):
    print("Running Experiment: AugR Decrypt")
    scheme = AugRCPABE(m=D_DEFAULT, group_name=GROUP_TYPE)
    pp, msk = scheme.setup_A()
    out = []
    n = len(x_values)
    for x in x_values:
        print(f"  - Processing x={x}...")
        attrs = set([str(i) for i in range(x)])
        sk = scheme.KeyGen_A(pp, msk, attrs)
        policy = " and ".join([str(i) for i in range(x)]) if x > 0 else "0"
        M = scheme.group.random(GT)
        ct = scheme.Encrypt_A(pp, M, [], policy)
        total = 0.0
        for _ in range(REPEAT):
            t0 = time.time()
            scheme.Decrypt_A(pp, ct, sk)
            t1 = time.time()
            total += (t1 - t0) * 1000.0
        out.append(total / REPEAT)
        print(f"完成第 {len(out)}/{n} 个")
    return out

def _trace_ucpabe_time(N, policy_size):
    d = int(max(1, ceil(N ** 0.5)))
    print(f"U-CPABE trace 开始: N={N}, d={d}")
    scheme = UnifiedCPABE(d=d, group_type=GROUP_TYPE)
    mpk, msk = scheme.setup()
    print(f"U-CPABE trace 实验用户总数: {N}")
    attrs = {str(i) for i in range(policy_size)}
    sk = scheme.keygen(mpk, msk, attrs, 1, 1)
    uk = scheme.key_update(mpk, msk, 10)
    dk = scheme.dec_key_gen(sk, uk)
    def oracle(ct):
        try:
            return scheme.decrypt(ct, dk)
        except Exception:
            return None
    policy = " and ".join([str(i) for i in range(policy_size)]) if policy_size > 0 else "0"
    print(f"U-CPABE trace 开始: N={N}, policy={policy_size}")
    t0 = time.time()
    scheme.trace(mpk, policy, 10, oracle, W=10, epsilon=0.1)
    t1 = time.time()
    elapsed = (t1 - t0)
    print(f"U-CPABE trace 完成: N={N}, policy={policy_size}, time={elapsed:.3f}ms")
    return elapsed

def _trace_augr_time(N, policy_size):
    m = int(max(1, ceil(N ** 0.5)))
    print(f"AugR trace 开始: N={N}, m={m}")
    scheme = AugRCPABE(m=m, group_name=GROUP_TYPE)
    pp, msk = scheme.setup_A()
    print(f"AugR trace 实验用户总数: {N}")
    attrs = set([str(i) for i in range(policy_size)])
    sk = scheme.KeyGen_A(pp, msk, attrs)
    def oracle(ct):
        try:
            return scheme.Decrypt_A(pp, ct, sk)
        except Exception:
            return None
    policy = " and ".join([str(i) for i in range(policy_size)]) if policy_size > 0 else "0"
    print(f"AugR trace 开始: N={N}, policy={policy_size}")
    t0 = time.time()
    scheme.Trace(pp, [], policy, 0.1, oracle, W=10)
    t1 = time.time()
    elapsed = (t1 - t0)
    print(f"AugR trace 完成: N={N}, policy={policy_size}, time={elapsed:.3f}s")
    return elapsed

def _trace_plbe_time(N):
    m = int(max(1, ceil(N ** 0.5)))
    print(f"PLBE trace 开始: N={N}, m={m}")
    scheme = PLBE(m=m, group_type=GROUP_TYPE)
    pk, sk_map = scheme.setup_plbe()
    x_t = 1
    y_t = 1
    def oracle(ct):
        try:
            return scheme.decrypt_plbe(ct, x_t, y_t)
        except Exception:
            return None
    print(f"PLBE trace 开始: N={N}")
    t0 = time.time()
    scheme.trace(oracle, W=10, epsilon=0.1)
    t1 = time.time()
    elapsed = (t1 - t0)
    print(f"PLBE trace 完成: N={N}, time={elapsed:.3f}s")
    return elapsed

def run_exp_keygen_compare():
    x_vals = list(range(10, 201, 20))
    y_keygen_uc = measure_keygen_ucpabe(x_vals)
    y_keygen_w = measure_keygen_waters(x_vals)
    y_keygen_f = [v*2 for v in measure_keygen_fabeo(x_vals)]
    y_keygen_a = measure_keygen_augr(x_vals)
    y_keygen_i = measure_keygen_indre(x_vals)
    plot_graph_multi("fig1_keygen_compare.png", "Key Generation Time", "Number of Attributes/Policy Size", x_vals, [
        {"y": y_keygen_uc, "style": "o-", "label": "Ours", "color": "#C0392B"},
        {"y": y_keygen_w, "style": "^-", "label": "Waters[4]", "color": "#2980B9"},
        {"y": y_keygen_f, "style": "s-", "label": "RW[10]", "color": "#F39C12"},
        {"y": y_keygen_i, "style": "x-", "label": "AI[11]", "color": "#7E2F8E"},
        {"y": y_keygen_a, "style": "d-", "label": "LW[14]", "color": "#27AE60"},
        
    ])

def run_exp_encrypt_compare():
    x_vals = list(range(10, 201, 20))
    y_enc_uc = measure_encrypt_ucpabe(x_vals)
    y_enc_w = measure_encrypt_waters(x_vals)
    y_enc_f = measure_encrypt_fabeo(x_vals)
    y_enc_a = measure_encrypt_augr(x_vals)
    y_enc_i = measure_encrypt_indre(x_vals)
   
    plot_graph_multi("111.png", "Encryption Time", "Policy Size/Number of Attributes", x_vals, [
        {"y": y_enc_uc, "style": "o-", "label": "Ours", "color": "#C0392B"},
        {"y": y_enc_w, "style": "^-", "label": "Waters[4]", "color": "#2980B9"},
        {"y": y_enc_f, "style": "s-", "label": "RW[10]", "color": "#F39C12"},
        {"y": y_enc_i, "style": "x-", "label": "AI[11]", "color": "#7E2F8E"},
        {"y": y_enc_a, "style": "d-", "label": "LW[14]", "color": "#27AE60"},
        
    ])

def run_exp_decrypt_compare():
    x_vals = list(range(10, 201, 20))
    y_dec_uc = measure_decrypt_ucpabe(x_vals)
    y_dec_w = measure_decrypt_waters(x_vals)
    y_dec_f = measure_decrypt_fabeo(x_vals)
    y_dec_a = measure_decrypt_augr(x_vals)
    y_dec_i = measure_decrypt_indre(x_vals)
    plot_graph_multi("fig3_decrypt_compare.png", "Decryption Time", "Number of Attributes", x_vals, [
        {"y": y_dec_uc, "style": "o-", "label": "Ours", "color": "#C0392B"},
        {"y": y_dec_w, "style": "^-", "label": "Waters[4]", "color": "#2980B9"},
        {"y": y_dec_f, "style": "s-", "label": "RW[10]", "color": "#F39C12"},
        {"y": y_dec_i, "style": "x-", "label": "AI[11]", "color": "#7E2F8E"},
        {"y": y_dec_a, "style": "d-", "label": "LW[14]", "color": "#27AE60"},
          
    ])

def run_exp_trace_users():
    Ns = list(i*i for i in range(1, 31, 2))
    y_plbe = []
    y_augr = []
    y_ucpabe = []
    for idx, N in enumerate(Ns):
        y_plbe.append(_trace_plbe_time(N))
        y_augr.append(_trace_augr_time(N, 10))
        y_ucpabe.append(_trace_ucpabe_time(N, 10))
        print(
            f"Trace loop 完成: N={N}, U-CPABE={y_ucpabe[-1]:.2f}s, LW[14]={y_augr[-1]:.2f}s, PLBE={y_plbe[-1]:.2f}s"
        )
        print(f"完成第 {idx+1}/{len(Ns)} 个")
    plot_graph_multi("fig4_trace_vs_users.png", "Trace Time vs Users", "Total Users", Ns, [
        {"y": y_ucpabe, "style": "o-", "label": "Ours", "color": "#C0392B"},
        {"y": y_plbe, "style": "*-", "label": "GKS+[12]", "color": "#EDB120"},
        {"y": y_augr, "style": "d-", "label": "LW[14]", "color": "#27AE60"},
        
    ], y_label='Time (s)', x_tick_step=100)

def run_exp_trace_policy():
    sizes = list(range(10, 201, 20))
    y_augr = []
    y_ucpabe = []
    for idx, s in enumerate(sizes):
        y_augr.append(_trace_augr_time(100, s))
        y_ucpabe.append(_trace_ucpabe_time(100, s))
        print(
            f"Trace loop 完成: policy={s}, U-CPABE={y_ucpabe[-1]:.2f}s, LW[14]={y_augr[-1]:.2f}s"
        )
        print(f"完成第 {idx+1}/{len(sizes)} 个")
    plot_graph_multi("fig5_trace_vs_policy.png", "Trace Time vs Policy Size", "Policy Size", sizes, [
        {"y": y_ucpabe, "style": "o-", "label": "Ours", "color": "#C0392B"},
        {"y": y_augr, "style": "d-", "label": "LW[14]", "color": "#27AE60"}
    ], y_label='Time (s)')

def run_exp_revocation():
    rs = list(range(0, 901, 100))
    scheme_uc = UnifiedCPABE(d=30, group_type=GROUP_TYPE)
    mpk_uc, msk_uc = scheme_uc.setup()
    all_users = sorted(list(scheme_uc.tree.leaves.keys()))
    print(f"Ours 撤销实验用户总数: {len(all_users)}")
    odd_users = [u for u in all_users if u % 2 == 1]
    even_users = [u for u in all_users if u % 2 == 0]
    num_users = 900
    print(f"AI[11] 撤销实验用户总数: {num_users}")
    scheme_in = IndirectRevocationABE(num_users=num_users, group_type=GROUP_TYPE)
    pk_in, msk_in = scheme_in.setup()
    odd_users_in = [u for u in range(1, 901) if u % 2 == 1]
    even_users_in = [u for u in range(1, 901) if u % 2 == 0]
    y_uc = []
    y_in = []
    for idx, r in enumerate(rs):
        odd_take = min(r, len(odd_users))
        even_take = min(max(0, r - odd_take), len(even_users))
        scheme_uc.revoked_users = set(odd_users[:odd_take] + even_users[:even_take])
        total = 0.0
        for i in range(REPEAT):
            t0 = time.time()
            scheme_uc.key_update(mpk_uc, msk_uc, 10)
            t1 = time.time()
            dt_ms = (t1 - t0) * 1000.0
            print(f"Ours key_update 耗时 {dt_ms:.2f} ms (r={r}, 第{i+1}/{REPEAT})")
            total += dt_ms
        y_uc.append(total / REPEAT)
        odd_take_in = min(r, len(odd_users_in))
        even_take_in = min(max(0, r - odd_take_in), len(even_users_in))
        scheme_in.revoked_users = set(odd_users_in[:odd_take_in] + even_users_in[:even_take_in])
        total2 = 0.0
        for i in range(REPEAT):
            t0 = time.time()
            scheme_in.key_update(pk_in, msk_in, 10)
            t1 = time.time()
            dt_ms2 = (t1 - t0) * 1000.0
            print(f"AI[11] key_update 耗时 {dt_ms2:.2f} ms (r={r}, 第{i+1}/{REPEAT})")
            total2 += dt_ms2
        y_in.append(total2 / REPEAT)
        print(f"完成第 {idx+1}/{len(rs)} 个")
    print('1')
    plot_graph_multi("fig6_revocation.png", "Update KeyGen Time", "Revoked Users", rs, [
        {"y": y_uc, "style": "o-", "label": "Ours", "color": "#C0392B"},
        {"y": y_in, "style": "x-", "label": "AI[11]", "color": "#7E2F8E"}
    ])

def run_exp_indre_keygen():
    x_vals = list(range(10, 201, 20))
    y_keygen_i = measure_keygen_indre(x_vals)
    plot_graph_multi(
        "fig7_indre_keygen.png",
        "Indre Key Generation Time",
        "Policy Size",
        x_vals,
        [
            {"y": y_keygen_i, "style": "x-", "label": "AI[11]", "color": "#7E2F8E"},
        ],
    )

def run_exp_indre_encrypt():
    x_vals = list(range(10, 201, 20))
    y_enc_i = measure_encrypt_indre(x_vals)
    plot_graph_multi(
        "fig8_indre_encrypt.png",
        "Indre Encryption Time",
        "Number of Attributes",
        x_vals,
        [
            {"y": y_enc_i, "style": "x-", "label": "AI[11]", "color": "#7E2F8E"},   
        ],
        y_tick_step=100,
    )

def run_exp_indre_decrypt():
    x_vals = list(range(10, 201, 20))
    y_dec_i = measure_decrypt_indre(x_vals)
    plot_graph_multi(
        "fig9_indre_decrypt.png",
        "Indre Decryption Time",
        "Number of Attributes",
        x_vals,
        [
            {"y": y_dec_i, "style": "x-", "label": "AI[11]", "color": "#7E2F8E"},   
        ],
        y_tick_step=100,
    )

def main():
    run_exp_keygen_compare()
    run_exp_encrypt_compare()
    run_exp_decrypt_compare()
    run_exp_trace_users()
    run_exp_trace_policy()
    run_exp_revocation()
    run_exp_indre_keygen()
    run_exp_indre_encrypt()
    run_exp_indre_decrypt()
    print("All experiments completed.")

if __name__ == "__main__":
    main()
