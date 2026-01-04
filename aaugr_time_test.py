import time
from charm.toolbox.pairinggroup import GT
from aaugr_cpabe import AugRCPABE

REPEAT = 5

def measure_keygen(x_values):
    scheme = AugRCPABE(m=10, group_name='MNT224')
    pp, msk = scheme.setup_A()
    print("KeyGen_A Time")
    for x in x_values:
        attrs = {str(i) for i in range(x)}
        total = 0.0
        for _ in range(REPEAT):
            t0 = time.time()
            scheme.KeyGen_A(pp, msk, attrs)
            t1 = time.time()
            total += (t1 - t0) * 1000.0
        print(f"x={x}, avg_ms={total/REPEAT:.2f}")

def measure_encrypt(x_values):
    scheme = AugRCPABE(m=10, group_name='MNT224')
    pp, msk = scheme.setup_A()
    print("Encrypt_A Time")
    for x in x_values:
        policy = " and ".join([str(i) for i in range(x)]) if x > 0 else "0"
        M = scheme.group.random(GT)
        total = 0.0
        for _ in range(REPEAT):
            t0 = time.time()
            scheme.Encrypt_A(pp, M, set(), policy)
            t1 = time.time()
            total += (t1 - t0) * 1000.0
        print(f"x={x}, avg_ms={total/REPEAT:.2f}")

def measure_decrypt(x_values):
    scheme = AugRCPABE(m=10, group_name='MNT224')
    pp, msk = scheme.setup_A()
    print("Decrypt_A Time")
    for x in x_values:
        attrs = {str(i) for i in range(x)}
        sk = scheme.KeyGen_A(pp, msk, attrs)
        policy = " and ".join([str(i) for i in range(x)]) if x > 0 else "0"
        M = scheme.group.random(GT)
        ct = scheme.Encrypt_A(pp, M, set(), policy)
        total = 0.0
        for _ in range(REPEAT):
            t0 = time.time()
            scheme.Decrypt_A(pp, ct, sk)
            t1 = time.time()
            total += (t1 - t0) * 1000.0
        print(f"x={x}, avg_ms={total/REPEAT:.2f}")

def main():
    xs = list(range(10, 201, 10))
    measure_keygen(xs)
    measure_encrypt(xs)
    measure_decrypt(xs)

if __name__ == "__main__":
    main()
