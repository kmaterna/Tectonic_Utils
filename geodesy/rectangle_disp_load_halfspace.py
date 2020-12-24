#!/usr/bin/env python
"""
Computing the horizontal and vertical displacements (u, v, w) at (x, y, z)
on an elastic half-space from a rectangle of uniform pressure p.
Rectangle is centered around the origin, and has x-dimensions 2a and y-dimensions 2b.
Kathryn Materna, June 21, 2017
"""

import numpy as np 
import matplotlib.pyplot as plt


def r_ij(i,j,x,y,z,a,b):
	"""
	r10 and r20: This is the radius between the point xyz and 
			the points x=a,y=0,z=0 or x=-a,y=0,z=0 
	r01 and r02: This is the radius between the point xyz and 
			the points x=0,y=b,z=0 or x=0,y=-b,z=0 
	"""
	if i==1 and j==0:
		r=np.sqrt((a-x)**2 + y**2 + z**2);
	elif i==2 and j==0:
		r=np.sqrt((a+x)**2 + y**2 + z**2);
	elif i==0 and j==1:
		r=np.sqrt(x**2 + (b-y)**2 + z**2);
	elif i==0 and j==2:
		r=np.sqrt(x**2 + (b+y)**2 + z**2);
	return r;
def beta_ij(i,j,x,y,z,a,b):
	"""
	This is a projection of the radius into one of the central planes.
	beta10 and beta20 are (a-/+x)^2 + z^2
	beta01 and beta02 are (b-/+y)^2 + z^2
	"""
	if j==0:   # the 10 or 20 cases
		beta=np.sqrt(r_ij(i,j,x,y,z,a,b)**2 - y**2);
	if j>0:    # the 01 or 02 cases
		beta=np.sqrt(r_ij(i,j,x,y,z,a,b)**2 - x**2);
	return beta;
def phi_ij(i,j,x,y,z,a,b):
	"""
	phi10 and phi20 are fractional radius with y on top.
	phi01 and phi02 are fractional radius with x on top. 
	"""
	tol=1e-9;
	if (r_ij(i,j,x,y,z,a,b)+beta_ij(i,j,x,y,z,a,b))<tol:
		phi=1;   # prevent dividing by zero; replace with the largest possible number. 
	else:
		if j==0:  # the 10 or 20 cases
			phi=y/(r_ij(i,j,x,y,z,a,b)+beta_ij(i,j,x,y,z,a,b));
		if j>0:  # the 01 or 02 cases
			phi=x/(r_ij(i,j,x,y,z,a,b)+beta_ij(i,j,x,y,z,a,b));
	return phi;

def J_j(j,x,y,z,a,b):
	tol=1e-9;
	if j==1:
		prefactor=np.abs(a-x);
	if j==2:
		prefactor=np.abs(a+x);
	if abs(z)<tol:  # limit issues with solving for surface displacements
		if abs(z+r_ij(j,0,x,y,z,a,b)) < tol:
			term1=0;
		else:
			term1=y*(np.log(z+r_ij(j,0,x,y,z,a,b))-1.0);
		term2=0;
		term3=2*prefactor*np.arctan(phi_ij(j,0,x,y,z,a,b));
	else:	
		term1=y*(np.log(z+r_ij(j,0,x,y,z,a,b))-1.0);
		term2=z*np.log((1.0+phi_ij(j,0,x,y,z,a,b))/(1.0-phi_ij(j,0,x,y,z,a,b)));
		term3=2*prefactor*np.arctan(prefactor*phi_ij(j,0,x,y,z,a,b)/(z+beta_ij(j,0,x,y,z,a,b)));
	J=term1+term2+term3;
	return J;

def K_j(j,x,y,z,a,b):
	tol=1e-9;
	if j==1:
		prefactor=np.abs(b-y);
	if j==2:
		prefactor=np.abs(b+y);
	if abs(z)<tol:  # limit issues with solving for surface displacements. 
		if abs(z+r_ij(0,j,x,y,z,a,b))<tol:
			term1=0;
		else:
			term1 = x*(np.log(z+r_ij(0,j,x,y,z,a,b))-1.0);
		term2=0;
		term3 = 2*prefactor*np.arctan(phi_ij(0,j,x,y,z,a,b));
	else:
		term1 = x*(np.log(z+r_ij(0,j,x,y,z,a,b))-1.0);
		term2 = z*np.log((1+phi_ij(0,j,x,y,z,a,b))/(1-phi_ij(0,j,x,y,z,a,b)));
		term3 = 2*prefactor*np.arctan(prefactor*phi_ij(0,j,x,y,z,a,b)/(z+beta_ij(0,j,x,y,z,a,b)));
	K = term1+term2+term3;
	return K;

def L_j(j,x,y,z,a,b):
	tol=1e-9;
	if j==1:
		prefactor=a-x;
	if j==2:
		prefactor=-a-x;
	if abs(z)<tol:
		if abs(prefactor+r_ij(j,0,x,y,z,a,b))<tol:
			term1=0;
		else:
			term1=y*(np.log(prefactor+r_ij(j,0,x,y,z,a,b))-1);
		if abs(prefactor)<tol:
			term2=0;
		else:
			term2=prefactor*np.log((1+phi_ij(j,0,x,y,z,a,b))/(1-phi_ij(j,0,x,y,z,a,b)));
		term3=0;
	else:
		term1=y*(np.log(prefactor+r_ij(j,0,x,y,z,a,b))-1);
		term2=prefactor*np.log((1+phi_ij(j,0,x,y,z,a,b))/(1-phi_ij(j,0,x,y,z,a,b)));
		term3=2*z*np.arctan((z*phi_ij(j,0,x,y,z,a,b))/(prefactor+beta_ij(j,0,x,y,z,a,b)));
	L = term1+term2+term3;
	return L;


def u(x, y, z, lame, mu, p, a, b):
	"""
	Displacements in the x direction. For the corners of the box, we set the exploding terms manually to 0.
	"""
	tol=1e-9;
	deltay_top = -b-y;
	deltay_bottom = b-y;
	if abs(deltay_top)<tol or abs(deltay_bottom)<tol or abs(-a-x)<tol or abs(a-x)<tol:
		if abs(z)<tol:
			u_top=(-p/4/np.pi)*((J_j(2,x,deltay_top,z,a,b)-J_j(1,x,deltay_top,z,a,b))*(1/(lame+mu)));
			u_bottom=(-p/4/np.pi)*((J_j(2,x,deltay_bottom,z,a,b)-J_j(1,x,deltay_bottom,z,a,b))*(1/(lame+mu)));
	else:
		u_top=(-p/4/np.pi)*((J_j(2,x,deltay_top,z,a,b)-J_j(1,x,deltay_top,z,a,b))*(1/(lame+mu))+(z/mu)*(np.log((deltay_top+r_ij(2,0,x,deltay_top,z,a,b))/(deltay_top+r_ij(1,0,x,deltay_top,z,a,b)))));
		u_bottom=(-p/4/np.pi)*((J_j(2,x,deltay_bottom,z,a,b)-J_j(1,x,deltay_bottom,z,a,b))*(1/(lame+mu))+(z/mu)*(np.log((deltay_bottom+r_ij(2,0,x,deltay_bottom,z,a,b))/(deltay_bottom+r_ij(1,0,x,deltay_bottom,z,a,b)))));
	return -1*(u_top - u_bottom);

def v(x, y, z, lame, mu, p, a, b):
	"""
	Displacements in the x direction. For the corners of the box, we set the exploding terms manually to 0.
	"""
	tol=1e-9;
	deltax_top = -a-x;
	deltax_bottom = a-x;
	if abs(deltax_top)<tol or abs(deltax_bottom)<tol or abs(-b-y)<tol or abs(b-y)<tol:
		v_top=(-p/4/np.pi)*((K_j(2,deltax_top,y,z,a,b)-K_j(1,deltax_top,y,z,a,b))*(1/(lame+mu)));
		v_bottom=(-p/4/np.pi)*((K_j(2,deltax_bottom,y,z,a,b)-K_j(1,deltax_bottom,y,z,a,b))*(1/(lame+mu)));
	else:
		v_top=(-p/4/np.pi)*((K_j(2,deltax_top,y,z,a,b)-K_j(1,deltax_top,y,z,a,b))*(1/(lame+mu))+(z/mu)*(np.log((deltax_top+r_ij(0,2,deltax_top,y,z,a,b))/(deltax_top+r_ij(0,1,deltax_top,y,z,a,b)))));
		v_bottom=(-p/4/np.pi)*((K_j(2,deltax_bottom,y,z,a,b)-K_j(1,deltax_bottom,y,z,a,b))*(1/(lame+mu))+(z/mu)*(np.log((deltax_bottom+r_ij(0,2,deltax_bottom,y,z,a,b))/(deltax_bottom+r_ij(0,1,deltax_bottom,y,z,a,b)))));
	return -1*(v_top - v_bottom);

def w(x,y,z,lame,mu,p,a,b):
	tol=1e-9;
	deltay_top = -b-y;
	deltay_bottom = b-y;
	if abs(z)<tol:
		w_top=(p/4/np.pi/mu)*((lame+2*mu)/(lame+mu))*(L_j(1,x,deltay_top,z,a,b)-L_j(2,x,deltay_top,z,a,b));
		w_bottom=(p/4/np.pi/mu)*((lame+2*mu)/(lame+mu))*(L_j(1,x,deltay_bottom,z,a,b)-L_j(2,x,deltay_bottom,z,a,b));
	else:	
		w_top=(p/4/np.pi/mu)*(((lame+2*mu)/(lame+mu))*(L_j(1,x,deltay_top,z,a,b)-L_j(2,x,deltay_top,z,a,b)) + z*(np.arctan(((a-x)*deltay_top)/(z*r_ij(1,0,x,deltay_top,z,a,b)))+np.arctan(((a+x)*deltay_top)/(z*r_ij(2,0,x,deltay_top,z,a,b)))))
		w_bottom=(p/4/np.pi/mu)*(((lame+2*mu)/(lame+mu))*(L_j(1,x,deltay_bottom,z,a,b)-L_j(2,x,deltay_bottom,z,a,b)) + z*(np.arctan(((a-x)*deltay_bottom)/(z*r_ij(1,0,x,deltay_bottom,z,a,b)))+np.arctan(((a+x)*deltay_bottom)/(z*r_ij(2,0,x,deltay_bottom,z,a,b)))))
	return (w_top - w_bottom);


if __name__=="__main__":
	a=500.0;
	b=1000.0;
	lame=24e9;
	mu=24e9; 
	depth_of_water = 100;
	p=1000*9.81*depth_of_water;  # 100m water depth

	x=np.array(np.arange(-2000.0,2000.0,1));
	u_displacement=[1000.0 * u(i, 0.0, 0.0, lame, mu, p, a, b) for i in x];
	v_displacement=[1000.0 * v(i, 0, 0.0, lame, mu, p, a, b) for i in x];
	w_displacement=[1000.0 * w(i, 0, 0.0, lame, mu, p, a, b) for i in x];
	# w_displacement=[1000.0 * w(500.0, 0, 0.0, lame, mu, p, a, b) ];

	# fig=plt.figure(figsize=(20,8));
	# f, axarr = plt.subplots(5, figsize=(20,8))
	fig=plt.figure(figsize=(20,8));
	ax1 = plt.subplot2grid((4, 2), (0, 0))
	ax1.plot(x,u_displacement);
	ax1.plot(x,v_displacement);
	ax1.text(-2000,0,'A',fontsize=18)
	ax1.text(2000,0,'A\'',fontsize=18)
	ax1.set_ylabel('H. Disp. (mm)')
	ax1.set_title('Elastic Displacement from '+str(depth_of_water)+'m EWT Load')

	ax2 = plt.subplot2grid((4, 2), (1, 0))
	ax2.plot(x,w_displacement);
	ax2.set_ylabel('V. Disp. (mm)')

	y=np.array(np.arange(-2000.0,2000.0,1));
	u_displacement=[1000.0 * u(0.0, i, 0.0, lame, mu, p, a, b) for i in y];
	v_displacement=[1000.0 * v(0.0, i, 0.0, lame, mu, p, a, b) for i in y];
	w_displacement=[1000.0 * w(0.0, i, 0.0, lame, mu, p, a, b) for i in y];

	ax3 = plt.subplot2grid((4, 2), (2, 0))
	ax3.plot(y,u_displacement);
	ax3.plot(y,v_displacement);
	ax3.text(-2000,0,'B',fontsize=18)
	ax3.text(2000,0,'B\'',fontsize=18)
	ax3.set_ylabel('H. Disp. (mm)')

	ax4 = plt.subplot2grid((4, 2), (3, 0))
	ax4.plot(y,w_displacement);
	ax4.set_ylabel('V. Disp. (mm)')
	ax4.set_xlabel('Distance (m)')

	
	ax5 = plt.subplot2grid((4, 2), (0, 1), rowspan=4)
	ax5.grid('on');
	ax5.plot([a, a, -a, -a, a], [b, -b, -b, b, b],'k');
	ax5.set_xlim([-2000, 2000])
	ax5.set_ylim([-2000, 2000])
	ax5.plot([-2000,2000],[0,0],linewidth=2,color='black');
	ax5.plot([0,0],[-2000,2000],linewidth=2,color='black');
	ax5.text(2000,0,'A\'',fontsize=24);
	ax5.text(-2000,0,'A',fontsize=24);
	ax5.text(0,-2000,'B',fontsize=24);
	ax5.text(0,2000,'B\'',fontsize=24);

	plt.savefig('u_disp.eps');


# Right now I have a little bug relating to the phi in the corners of the rectangle. 
