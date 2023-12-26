import numpy as np

def inner_surface(x, y, z, r, rn, chn, delz):
    """
    Extracts inner surface data and writes it to a CSV file.

    Parameters:
    - x, y, z: Arrays representing coordinates.
    - r: Array representing radius.
    - rn: Array representing something related to radius (e.g., not specified in the original code).
    - chn: Array representing chain information.
    - delz: Step size for the z-axis.

    Returns:
    - All_data: List containing arrays for x, y, z, r, rn, and chn of the inner surface.
    """

    # Open the CSV file for writing
    with open("1C8F_T1.csv", "w") as f:
        z_min = min(z)
        z_max = max(z)
        Zr = np.arange(np.floor(z_min), np.ceil(z_max), delz)

        x_main, y_main, z_main, r_main, rn_main, chn_main = [], [], [], [], [], []

        for i in range(len(Zr) - 1):
            x_zr, y_zr, z_zr, r_zr, rn_zr, chn_zr = [], [], [], [], [], []
            first_rad = 0

            for j in range(len(z)):
                if Zr[i] <= z[j] < Zr[i + 1]:
                    x_zr.append(x[j])
                    y_zr.append(y[j])
                    z_zr.append(z[j])
                    r_zr.append(r[j])
                    rn_zr.append(rn[j])
                    chn_zr.append(chn[j])

            rad_store = [np.sqrt(x**2 + y**2) for x, y in zip(x_zr, y_zr)]
            first_rad = min(rad_store)

            for j in range(len(x_zr)):
                if first_rad**2 <= x_zr[j]**2 + y_zr[j]**2 < (10 + first_rad)**2:
                    x_main.append(x_zr[j])
                    y_main.append(y_zr[j])
                    z_main.append(z_zr[j])
                    r_main.append(r_zr[j])
                    rn_main.append(rn_zr[j])
                    chn_main.append(chn_zr[j])
                    f.write(f'{r_zr[j]} \t {rn_zr[j]} \t {chn_zr[j]}\n')

        All_data = [np.array(x_main), np.array(y_main), np.array(z_main),
                    np.array(r_main), np.array(rn_main), np.array(chn_main)]

    return All_data


