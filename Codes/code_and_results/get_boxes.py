from mdtraj.formats import XTCTrajectoryFile


xtcname = 'water_traj/bulk_tip3p_all_part1.xtc'
traj_file = XTCTrajectoryFile(xtcname)
_, _ , _,box = traj_file.read()

file1 = open("boxes.txt", 'w')

for item in box:
    file1.write(f'{item[0][0]}, {item[1][1]}, {item[2][2]}\n')

file1.close()