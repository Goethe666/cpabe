
from charm.toolbox.pairinggroup import GT
from u_cpabe import UnifiedCPABE
import time


def test_basic_decrypt():
    t1_list=[]
    t2_list=[]
    t3_list=[]
    t4_list=[]
    t5_list=[]
    for i in range (10):

        scheme = UnifiedCPABE(d=2, group_type='MNT224')
        mpk, msk = scheme.setup()
        policy = "(EMPLOYEE and FINANCE) or (1 and 2) or 3"
        t = 2
        M = scheme.group.random(GT)
        start_time = time.time()
        sk = scheme.keygen(mpk, msk, {"EMPLOYEE", "FINANCE", 1, 2, 3}, 1, 1)
        end_time = time.time()
        t1=end_time - start_time
        t1_list.append(t1)  
        print("密钥生成时间:", t1)
        start_time = time.time()
        uk = scheme.key_update(mpk, msk, t)
        end_time = time.time()
        t2=end_time - start_time
        t2_list.append(t2)  
        print("密钥更新时间:", t2)
        start_time=time.time()
        dk = scheme.dec_key_gen(sk, uk)
        end_time=time.time()
        t3=end_time - start_time
        t3_list.append(t3)          
        print("解密密钥生成时间:", t3)
        start_time=time.time()
        ct = scheme.encrypt(mpk, M, policy, t)
        end_time=time.time()
        t4=end_time - start_time
        t4_list.append(t4)          
        print("加密时间:", t4)
        start_time=time.time()
        rec = scheme.decrypt(ct, dk)
        
        print("解密结果:", rec)
        end_time=time.time()
        t5=end_time - start_time
        t5_list.append(t5)          
        print("解密时间:", t5)
        print("基础解密:", "OK" if rec == M else "FAIL")
    s1=sum(t1_list)/10
    s2=sum(t2_list)/10
    s3=sum(t3_list)/10
    s4=sum(t4_list)/10
    s5=sum(t5_list)/10
    print("密钥生成时间:", s1)
    print("密钥更新时间:", s2)
    print("解密密钥生成时间:", s3)
    print("加密时间:", s4)
    print("解密时间:", s5)
def test_revocation():
    scheme = UnifiedCPABE(d=2, group_type='MNT224')
    mpk, msk = scheme.setup()
    policy = "(EMPLOYEE and FINANCE)or (1 and 2) or 3"
    t = 2
    target = (1, 2)
    sk_user = scheme.keygen(mpk, msk, {"EMPLOYEE", "FINANCE",1,2,3}, target[0], target[1])
    scheme.revoked_users.add(scheme._user_id(*target))
    uk = scheme.key_update(mpk, msk, t)
    dk = scheme.dec_key_gen(sk_user, uk)
    if dk is None:
        print("撤销机制: OK (DecKeyGen 返回 None)")
    else:
        M = scheme.group.random(GT)
        ct = scheme.encrypt(mpk, M, policy, t)
        rec = scheme.decrypt(ct, dk)
        print("撤销解密:", "BLOCKED" if rec != M else "FAIL")

def test_trace(traitor_pos=(2, 1), epsilon=0.4, W=30):
    scheme = UnifiedCPABE(d=2, group_type='MNT224')
    mpk, msk = scheme.setup()
    policy = "(EMPLOYEE and FINANCE)"
    t = 2
    sk_traitor = scheme.keygen(mpk, msk, {"EMPLOYEE", "FINANCE"}, traitor_pos[0], traitor_pos[1])
    uid = scheme._user_id(*traitor_pos)
    if uid in scheme.revoked_users:
        scheme.revoked_users.remove(uid)
    uk = scheme.key_update(mpk, msk, t)
    dk_traitor = scheme.dec_key_gen(sk_traitor, uk)
    if dk_traitor is None:
        print("追踪: 密钥生成失败")
        return
    def decode_oracle(ct):
        try:
            return scheme.decrypt(ct, dk_traitor)
        except Exception:
            return None
    detected = scheme.trace(mpk, policy, t, decode_oracle, W=W, epsilon=epsilon)
    print("追踪结果:", detected)
    print("追踪校验:", "OK" if traitor_pos in detected else "FAIL")

def main():
    test_basic_decrypt()
    test_revocation()
    test_trace((2, 1), 0.4, 30)

if __name__ == "__main__":
    main()
