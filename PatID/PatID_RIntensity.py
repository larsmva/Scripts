


def mric_RI(save , *args):

	import numpy as np
	import nibabel as nb
	
	target = args[0][0]
	
	vox2ras = nb.load(target).get_header().get_vox2ras()

	idata  = nb.load(target).get_data()
	init = idata.reshape(1,-1)[0].astype("float32")

	for j in range(1,len(args[0])):
		i = args[0][j]
		
		temp=nb.load(i).get_data().reshape(1,-1)[0].astype("float32")
	    
		rI2 = temp/init
		rI2[np.isnan(rI2)] = 1 #
		temp3 = rI2.reshape(idata.shape)
		img2 = nb.MGHImage(temp3.astype("float32"),vox2ras)
	
		print save+"/"+str(i).split("/")[-1]
		
		nb.save(img2,save+"/"+str(i).split("/")[-1])
		
if __name__ =='__main__':
	import sys 

	mric_RI(sys.argv[1] , sys.argv[2::])
