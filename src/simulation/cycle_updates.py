import src.model.firing_rate as firing
from src.utils.helpers import get_c

def get_alpha(n_forks, chi, volume, regime):
    """
        returns the fraction of active DnaA
    """
    chi=n_forks*chi
    if regime=="constant":
        alpha=(volume-chi)/volume
    elif regime=="linear":
        alpha=volume/(volume+chi)
    else:
        print("unknown regime in get_alpha_opt")
        return 0
    return max((alpha, 1e-10))

def update_volume(volume, dt, cfg):
    """
        Exponentially updates cell volume based on the growth rate.
    """
    return volume*(1. + cfg.model.GROWTH_RATE*dt)

def update_n_titration(n_tot, n_forks, dt, cfg):
    """
        Updates the total number of titration sites over time.

        Assumes that titration sites are uniformly distributed along the chromosome.
        As replication progresses, new sites are added proportionally to the number of forks.

    """
    return n_tot + n_forks*cfg.model.SITES/(2*cfg.model.REP_TIME)*dt

def make_step(n_forks, n_tot, volume, a_atp, a_adp, time, dt, y, chi, step, cfg):
    """
        This function updates the cell volume, titration site count, DnaA activation state,
        and computes the firing rate. 
    """
    volume = update_volume(volume, dt, cfg)
    n_tot = update_n_titration(n_tot, n_forks, dt, cfg=cfg)
    c_tot = n_tot/volume
    alpha = get_alpha(n_forks, chi, volume, regime=cfg.model.REGIME)
    a_atp, a_adp = alpha*cfg.model.DNAA_CONCENTRATION, (1.-alpha)*cfg.model.DNAA_CONCENTRATION  
    c = get_c(cfg.model.DNAA_CONCENTRATION, cfg.model.K, c_tot)
    c_atp, c_adp = alpha*c, (1.-alpha)*c
    if step:
        #if a_atp>c_tot-c_adp:
        #if a_atp-c_atp>10.:
        if cfg.model.REGIME=="linear":
            if volume>(n_tot+10.*chi*n_forks)/(cfg.model.DNAA_CONCENTRATION-10.):
                f_rate=1e6
            else:
                f_rate=0.
        else: 
            if a_atp-c_atp>10.:
                f_rate=1e6
            else: 
                f_rate=0
    else:
        f_rate = firing.fr(a_atp, a_adp, c_atp, c_adp, y, 
                        kori=cfg.model.K_OPEN, 
                        ori_sites=cfg.model.ORIGIN_SITES, 
                        epsilon_cost=cfg.model.E_COST, 
                        k_max=cfg.model.FIRING_MAX) 
    time += dt
    return time, a_atp, a_adp, c_atp, c_adp, volume, n_tot, f_rate

def main():
    return 0

if __name__=="__main__":
    main()