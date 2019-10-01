import wells_and_coppersmith as wc

mag = 6.8
RLD = wc.RLD_from_M(mag, "N");
print("RLD = %f km" % RLD)

RW = wc.RW_from_M(mag, "N");
print("RW = %f km" % RW)



