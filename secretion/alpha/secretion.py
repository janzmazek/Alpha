from .montefusco_cy import montefusco_euler as secretion
import numpy as np
import multiprocessing

MAX_AUTOCRINE = None

def autocrine_glucagon(normalized):
    return (1 + 175*normalized**8/(2.8**8+normalized**8))*normalized

def montefusco(alpha, g, cAMP, AA, return_dict, first_iteration=True):
    print(g)
    g_L = lambda g_K_ATP: 0.213

    # Defining initial values
    t = 100000
    dt = 0.1 if first_iteration else 0.05
    x_0 = np.array([-5.14063726e+01,  1.05155470e-01,  7.76130678e-01,
                    3.36181980e-06, 7.76130678e-01,  3.68483231e-01,
                    3.97062254e-01, 4.67787961e-03, 2.09233487e-01,
                    2.81686443e-01,  3.45055157e-01, 1.54643301e-01,
                    2.37897318e-01,  2.94083254e-01,  9.85403353e+01]
                )
    
    gk = alpha.g_K_ATP(g)
    voltage1, currents1 = secretion(x_0, gk, g_L(gk),
                                    cAMP, AA, t, dt
                                    )
    GS = np.mean(currents1[int(4*t/5):, -1])
    V_min = np.min(voltage1[int(4*t/5):, 0])
    V_max = np.max(voltage1[int(4*t/5):, 0])

    if g == 0 and AA == 0 and first_iteration:
        return_dict["max_autocrine"] = GS

    return_dict[g] = (GS, V_min, V_max)


def run_montefusco_parallel(alpha, glucose, amino_acids):
    """
    parameters: 
    - alpha object
    - glucose (np.arange(0, MAX_GLUCOSE+DIFF, DIFF))
    - amino_acids ([0, 0.5, 1])
    """

    for AA in amino_acids:
        print(f"AA:{AA}")
        manager = multiprocessing.Manager()
        return_dict = manager.dict()

        sAC, tmAC = [], []

        pool = multiprocessing.Pool(12)

        # Start multiprocesses (iteration 1)
        for g in glucose:
            cAMP = 0.75*alpha.cAMP_sAC_interpolation(g)
            pool.apply_async(montefusco, args=(alpha, g, cAMP, AA, return_dict))

        # Wait for multiprocesses to finish (iteration 1)
        pool.close()
        pool.join()
        pool = multiprocessing.Pool(12)

        # Normalize GS of first iteration
        if "max_autocrine" in return_dict:
            MAX_AUTOCRINE = return_dict["max_autocrine"]
        GS1 = [return_dict[g][0]/MAX_AUTOCRINE for g in glucose]

        # Start multiprocesses (iteration 2)
        for i, g in enumerate(glucose):
            cAMP_sAC = alpha.cAMP_sAC_interpolation(g)
            sAC.append(cAMP_sAC)
            cAMP_tmAC = autocrine_glucagon(GS1[i])
            tmAC.append(cAMP_tmAC)
            cAMP = 0.75*cAMP_sAC + 0.25*cAMP_tmAC
            fAA = AA

            pool.apply_async(montefusco, args=(alpha, g, cAMP, fAA, return_dict, False))

        # Wait for multiprocesses to finish (iteration 2)
        pool.close()
        pool.join()

        GS2 = [return_dict[g][0] for g in glucose]

        result = [[g, GS1[i], GS2[i], sAC[i], tmAC[i]] for i, g in enumerate(glucose)]
        yield (AA, np.array(result))
