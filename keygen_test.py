import time
from u_cpabe import UnifiedCPABE
from aaugr_cpabe import AugRCPABE

REPEAT = 5
D_DEFAULT = 4


def test_ucpabe_keygen(x_values, group_type='BN254'):  
    print("U-CPABE KeyGen Test")
    for x in x_values:
        scheme = UnifiedCPABE(d=D_DEFAULT, group_type=group_type)
        mpk, msk = scheme.setup()
        attrs = {str(i) for i in range(x)}
        total = 0.0
        for _ in range(REPEAT):
            t0 = time.time()
            scheme.keygen(mpk, msk, attrs, 1, 1)
            t1 = time.time()
            total += (t1 - t0) * 1000.0
        print(f"x={x}, avg_ms={total/REPEAT:.2f}")

def test_augr_keygen(x_values):
    print("AugR-CPABE KeyGen Test")
    scheme = AugRCPABE(m=D_DEFAULT, group_name='BN254')
    pp, msk = scheme.setup_A()
    for x in x_values:
        attrs = {str(i) for i in range(x)}
        total = 0.0
        for _ in range(REPEAT):
            t0 = time.time()
            scheme.KeyGen_A(pp, msk, attrs)
            t1 = time.time()
            total += (t1 - t0) * 1000.0
        print(f"x={x}, avg_ms={total/REPEAT:.2f}")

def main():
    xs = list(range(10, 201, 20))
    test_ucpabe_keygen(xs)
    test_augr_keygen(xs)

if __name__ == "__main__":
    main()
