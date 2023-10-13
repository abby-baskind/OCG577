


def H_poly2(H, AT, DIC,K1,K2,KB,BT):
    value=(H*H*H*AT)+H*H*(AT*KB+AT*K1-KB*BT-DIC*K1)+H*(AT*KB*K1+AT*K1*K2-KB*BT*K1-K1*KB*DIC-2*K1*K2*DIC) \
        + (AT*KB*K1*K2-KB*K1*K2*BT - 2*K1*K2*KB*DIC)
    return value