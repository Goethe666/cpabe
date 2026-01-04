import time
import matplotlib.pyplot as plt

from u_cpabe import UnifiedCPABE
from charm.toolbox.pairinggroup import GT


REPEAT = 20 
D_DEFAULT = 2 
D_REVOCATION =16
GROUP_TYPE = 'MNT224'

def run_experiment(name, x_values, exp_type):
    print(f"Running Experiment: {name}")
    times = []
    for x in x_values:
        print(f"  - Processing x={x}...")
        total_ms = 0
        if exp_type == "KeyGen":
            scheme = UnifiedCPABE(d=D_DEFAULT, group_type=GROUP_TYPE)
            mpk, msk = scheme.setup()
            attrs = {str(i) for i in range(x)}
            for _ in range(REPEAT):
                t_start = time.time()
                scheme.keygen(mpk, msk, attrs, 1, 1)
                t_end = time.time() 
                total_ms += (t_end - t_start) * 1000.0
        elif exp_type == "KeyUpdate":
            scheme = UnifiedCPABE(d=D_REVOCATION, group_type=GROUP_TYPE)
            mpk, msk = scheme.setup()
            all_users = sorted(list(scheme.tree.leaves.keys()))
            odd_users = [u for u in all_users if u % 2 == 1]
            even_users = [u for u in all_users if u % 2 == 0]
            odd_take = min(x, len(odd_users))
            even_take = min(max(0, x - odd_take), len(even_users))
            scheme.revoked_users = set(odd_users[:odd_take] + even_users[:even_take])
            cover_nodes = scheme.tree.cover_algorithm(scheme.revoked_users)
            print("Minimal cover set (node_ids):", [n.node_id for n in cover_nodes])
            t = 10
            for _ in range(REPEAT):
                t_start = time.time()
                scheme.key_update(mpk, msk, t)
                t_end = time.time()
                total_ms += (t_end - t_start) * 1000.0
        elif exp_type == "Encrypt":
            scheme = UnifiedCPABE(d=D_DEFAULT, group_type=GROUP_TYPE)
            mpk, msk = scheme.setup()
            policy_str = " and ".join([str(i) for i in range(x)]) if x != 1 else "0"
            M = scheme.group.random(GT)
            t = 1
            for _ in range(REPEAT):
                t_start = time.time()
                scheme.encrypt(mpk, M, policy_str, t)
                t_end = time.time()
                total_ms += (t_end - t_start) * 1000.0
        elif exp_type == "Decrypt_KeyAttrs":
            scheme = UnifiedCPABE(d=D_DEFAULT, group_type=GROUP_TYPE)
            mpk, msk = scheme.setup()
            policy_size = x
            policy_str = " and ".join([str(i) for i in range(policy_size)])
            print(f"Policy: {policy_str}")
            attrs = {str(i) for i in range(x)}
            print(f"Attributes: {attrs}")
            sk = scheme.keygen(mpk, msk, attrs, 1, 1)
            M = scheme.group.random(GT)
            t = 10
            ct = scheme.encrypt(mpk, M, policy_str, t)
            uk = scheme.key_update(mpk, msk, t)
            dk = scheme.dec_key_gen(sk, uk)
            for _ in range(REPEAT):
                t_start = time.time()
                scheme.decrypt(ct, dk)
                t_end = time.time()

                total_ms += (t_end - t_start) * 1000.0
        avg_time = total_ms / REPEAT
        times.append(avg_time)
    return times

def plot_graph(filename, title, x_label, x_values, y_values, marker='ro-', label=None, y_top=None):
    plt.figure(figsize=(8, 6))
    plt.plot(x_values, y_values, marker, linewidth=2, markersize=8, label=label)
    plt.title(title, fontsize=14)
    plt.xlabel(x_label, fontsize=12)
    plt.ylabel('Time (ms)', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xlim(left=0)
    if y_top is not None:
        plt.ylim(bottom=0, top=y_top)
    else:
        plt.ylim(bottom=0)
    if len(x_values) >= 2:
        step = x_values[1] - x_values[0]
        if step > 0:
            xmax = max(x_values)
            plt.xticks(list(range(0, int(xmax) + int(step), int(step))))
    if label:
        plt.legend(shadow=True, loc='upper left', fontsize=10)
    plt.tight_layout()
    plt.savefig(filename)
    print(f"Saved {filename}")
    plt.close()

 

def main():
    x1 = list(range(10, 201, 10))
    y1 = run_experiment("KeyGen", x1, "KeyGen")
    plot_graph("fig1_KeyGen.png", "Key Generation Time", 
               "Number of Attributes", x1, y1, 'ro-', "U-CPABE")

    x2 = list(range(0, 201, 10))
    y2 = run_experiment("KeyUpdate", x2, "KeyUpdate")
    plot_graph("fig2_keyupdate.png", "Key Update Time", 
               "Revocation List Size", x2, y2, 'gx-', "U-CPABE")

    x3 = list(range(10, 201, 10))
    y3 = run_experiment("Encrypt", x3, "Encrypt")
    plot_graph("fig3_encrypt.png", "Encryption Time", 
               "Policy Size", x3, y3, 'bs-', "U-CPABE") 

    x4 = list(range(10, 201, 10))
    y4 = run_experiment("Decrypt_KeyAttrs", x4, "Decrypt_KeyAttrs")
    plot_graph("fig4_decrypt_key.png", "Decryption Time", 
               "Number of Attributes", x4, y4, 'md-', "U-CPABE")

    

    

    print("All experiments completed.")

if __name__ == "__main__":
    main()
