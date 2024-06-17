""" This code will calculate read counts for unmapped vs mapped. Please note files must be .fq and within a directory labeled bbmap"""
import subprocess
out_unmapped = subprocess.check_output("cat bbmap/*unmapped.fq | wc -l",shell=True)
out_mapped = subprocess.check_output("cat bbmap/*mapped.fq | wc -l",shell=True)
unmapp = int(out_unmapped)/4
mapp = int(out_mapped)/4
tot = unmapp + mapp

unmap_p = unmapp * 100 / tot
map_p = mapp * 100 /tot

print("The total counts is {}".format(tot))
print("The map counts is {}".format(mapp))
print("The unmapp counts is {}".format(unmapp))
print("The unmap percentage is {}".format(unmap_p))
print("The map percentage  is {}".format(map_p))
