from charm.toolbox.pairinggroup import PairingGroup, GT
from augr_cpabe import AugRCPABE
from charm.schemes.abenc.abenc_waters11 import CPabe_Waters11


def test_basic():
    scheme = AugRCPABE(m=3, group_name='SS512')
    pp, msk = scheme.setup_A()
    sk = scheme.KeyGen_A(pp, msk, {"EMPLOYEE", "FINANCE"})
    M = scheme.group.random(GT)
    ct = scheme.Encrypt_A(pp, M, set(), "(EMPLOYEE and FINANCE)")
    rec = scheme.Decrypt_A(pp, ct, sk)
    print("basic:", "OK" if rec == M else "FAIL")


def test_revocation():
    scheme = AugRCPABE(m=3, group_name='SS512')
    pp, msk = scheme.setup_A()
    sk = scheme.KeyGen_A(pp, msk, {"EMPLOYEE", "FINANCE"})
    M = scheme.group.random(GT)
    ct = scheme.Encrypt_A(pp, M, {(sk['i'], sk['j'])}, "(EMPLOYEE and FINANCE)")
    rec = scheme.Decrypt_A(pp, ct, sk)
    print("revocation:", "BLOCKED" if rec != M else "FAIL")


def test_mismatch():
    scheme = AugRCPABE(m=3, group_name='SS512')
    pp, msk = scheme.setup_A()
    sk = scheme.KeyGen_A(pp, msk, {"EMPLOYEE"})
    M = scheme.group.random(GT)
    ct = scheme.Encrypt_A(pp, M, set(), "(EMPLOYEE and FINANCE)")
    rec = scheme.Decrypt_A(pp, ct, sk)
    print("mismatch:", "BLOCKED" if rec is None else "FAIL")


def main():
    test_basic()
    test_revocation()
    test_mismatch()
    test_waters11()


def test_trace():
    scheme = AugRCPABE(m=3, group_name='SS512')
    pp, msk = scheme.setup_A()
    policy = "(EMPLOYEE and FINANCE)"
    R_D = set()

    traitor_k = 5
    sk_traitor = None
    for _ in range(traitor_k):
        sk_traitor = scheme.KeyGen_A(pp, msk, {"EMPLOYEE", "FINANCE"})

    def decoder_oracle(ct):
        return scheme.Decrypt_A(pp, ct, sk_traitor)

    epsilon = 0.1
    traitors_idx = scheme.Trace(pp, R_D, policy, epsilon, decoder_oracle, W=20)
    print("trace-idx: expected", traitor_k, "found", traitors_idx)


def test_waters11():
    group = PairingGroup('SS512')
    cpabe = CPabe_Waters11(group)
    msg = group.random(GT)
    attributes = ['EMPLOYEE', 'FINANCE']
    access_policy = '((EMPLOYEE) and (FINANCE))'
    mpk, msk = cpabe.setup()
    sk = cpabe.keygen(mpk, msk, attributes)
    ct = cpabe.encrypt(mpk, msg, access_policy)
    rec = cpabe.decrypt(mpk, sk, ct)
    print("waters11_builtin:", "OK" if rec == msg else "FAIL")


if __name__ == "__main__":
    main()
    test_trace()
