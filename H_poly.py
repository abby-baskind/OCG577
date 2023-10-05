


def H_poly(H, AT, CO2_star, k0,k1,k2,kb,BT):
    value=H**3*AT+H**2*(AT*kb-kb*BT-CO2_star*k1)+ H*(-kb*k1*CO2_star-2*k1*k2*CO2_star) \
        - 2*k1*k2*kb*CO2_star
    return value