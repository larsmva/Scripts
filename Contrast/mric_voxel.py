

def mric_brentq(k1,k2,G,theta,I,max_num=3000 , tol=1.0e-6): #TO DO improve, implement with c++ boos
	
	if G==0:
		return 0
 	else :
		import numpy as np
	
		import scipy.optimize as sp
		aux0 = np.cos(np.deg2rad(theta)) 
		aux1= 1.-aux0
		aux2 = 1.+aux0
		x = ( k1*aux1 + k2*aux2  + np.sqrt(k1**2*aux1**2 + 2*k1*k2*aux1*aux2 + k2**2*aux1**2))/(2*k2) # if /aux0
		c_max = np.log(x/G)/k1  	
		
		
	 	def RI(c) : 
			return (1.-aux0*G)*(1.-G*np.exp(-k1*c))*np.exp(-k2*c)/((1.-aux0*G*np.exp(-k1*c))*(1.-G))
		Im = RI(c_max) 

		if I>Im or I <1.15:
		   return 0 # I = 1.0
	

		def F(c) : 
			return (1.-G*np.exp(-k1*c))*np.exp(-k2*c)/(1.-aux0*G*np.exp(-k1*c))  -I*(1.-G)/(1.-aux0*G) 
	     

		c =  sp.brentq(F,0.0,c_max)
	
	
		return c*1000

def mric_solve(segmented, target, TR,TE, alpha, save , *args):

	import nibabel as nb
	import numpy as np
	
	#Flatten array

	fseg= nb.load(segmented).get_data().reshape(1,-1)[0]
	
	rI = np.zeros(len(fseg))


	vox2ras = nb.load(target).get_header().get_vox2ras()

	init  = nb.load(target).get_data().reshape(1,-1)[0].astype("float32")
	#init[init==0]=0.1 # To avoid warnings
	solver = np.vectorize(mric_brentq)

	G =  np.zeros(len(fseg))
	k1 =TR*3.2
	k2 =TE*3.9 
 	   
	G[fseg>0] =   np.exp(-TR/1.084)          #white  assume all white 
	G[fseg>999] =      np.exp(-TR/1.820)         #grey
        G[fseg==4]  = 0  # np.exp(-TR/20.000) 
	G[fseg==43] =  0 #   np.exp(-TR/20.000)      #ve 


	for j in range(len(args)):
		i = args[j]

		temp=nb.load(i).get_data().reshape(1,-1)[0].astype("float32")
	    
		rI[fseg>0] = temp[fseg>0]/init[fseg>0]
		

		
		

		rI[np.isnan(rI)] = 1 # can be included 
	

		

		temp2 = np.array(solver(k1,k2,G,alpha,rI)).reshape(256,256,256)
		img = nb.MGHImage(temp2.astype("float32"),vox2ras)
		nb.save(img,save+str(j)+".mgz")


		rI2 = temp/init
		rI2[np.isnan(rI2)] = 1 #
		temp3 = rI2.reshape(256,256,256)
		img2 = nb.MGHImage(temp3.astype("float32"),vox2ras)
		nb.save(img2,"testIR"+str(j)+".mgz")


def mric_RI(target, save , *args):

	import numpy as np
	import nibabel as nb

	vox2ras = nb.load(target).get_header().get_vox2ras()

	idata  = nb.load(target).get_data()
	init = idata.reshape(1,-1)[0].astype("float32")
	
	
	for j in range(len(args)):
		i = args[j]

		temp=nb.load(i).get_data().reshape(1,-1)[0].astype("float32")
	    
		rI2 = temp/init
		rI2[np.isnan(rI2)] = 1 #
		temp3 = rI2.reshape(idata.shape)
		img2 = nb.MGHImage(temp3.astype("float32"),vox2ras)
		print str(i).rsplit("/")[1]
		nb.save(img2,save+"RI"+str(i).rsplit("/")[1])
		
