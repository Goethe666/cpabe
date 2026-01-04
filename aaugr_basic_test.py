from charm.toolbox.pairinggroup import GT
from aaugr_cpabe import AugRCPABE

def test_basic():
    scheme = AugRCPABE(m=3, group_name='MNT224')
    pp, msk = scheme.setup_A()
    sk = scheme.KeyGen_A(pp, msk, {"EMPLOYEE", "FINANCE"})
    M = scheme.group.random(GT)
    ct = scheme.Encrypt_A(pp, M, set(), "(EMPLOYEE and FINANCE)")
    rec = scheme.Decrypt_A(pp, ct, sk)
    print("basic:", "OK" if rec == M else "FAIL")

def test_mismatch():
    scheme = AugRCPABE(m=3, group_name='MNT224')
    pp, msk = scheme.setup_A()
    sk = scheme.KeyGen_A(pp, msk, {"EMPLOYEE"})
    M = scheme.group.random(GT)
    ct = scheme.Encrypt_A(pp, M, set(), "(EMPLOYEE and FINANCE)")
    rec = scheme.Decrypt_A(pp, ct, sk)
    print("mismatch:", "BLOCKED" if rec is None else "FAIL")

def main():
    test_basic()
    test_mismatch()

if __name__ == "__main__":
    main()
