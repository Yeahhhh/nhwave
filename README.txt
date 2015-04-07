a 2015 GPU Hackathon code
	based on nhwave2.0.tar.gz
	port to GPU with OpenACC


reference
	https://www.olcf.ornl.gov/training-event/2015-gpu-hackathons/	2015 GPU Hackathons
	http://on-demand.gputechconf.com/gtc/2015/video/S5515.html	A summary of the first GPU hackathon



===================================================

Original Statements
https://sites.google.com/site/gangfma/nhwave


NHWAVE

	NHWAVE is a three-dimensional shock-capturing Non-Hydrostatic WAVE
	model developed by Ma et al. (2012), which solves the incompressible
	Navier-Stokes equations in terrain and surface-following sigma coordinates.
	The model predicts instantaneous surface elevation and 3D flow field, and is
	capable of resolving coastal wave processes (shoaling, refraction,
	diffraction, breaking etc.) as well as tsunami wave generation by submarine
	mass failure. The governing equations are discretized by a shock-capturing
	Godunov-type numerical scheme. A nonlinear Strong Stability-Preserving (SSP)
	Runge-Kutta scheme is adopted for adaptive time stepping with second-order
	temporal accuracy. The model is fully parallelized using Message Passing
	Interface (MPI) with non-blocking communication. The poisson equation is
	solved by the high performance preconditioner HYPRE software library
	(http://acts.nersc.gov/hypre/). The details of the model are referred to Ma et
	al. (2012). 

Major Features:
	Model solves 3D Navier-Stokes equations in surface and terrain-following sigma coordinate;
	Godunov-type shock-capturing TVD scheme was employed;
	Adaptive time stepping was adopted using second-order nonlinear Strong
	Stability-Preserving (SSP) Runge-Kutta scheme;
	Adaptive grid refinement (AMR) technique was implemented;
	The model is capable of simulating tsunami wave generation by submarine landslide (bottom movement);
	A cohesive/non-cohesive sediment transport module was implemented taking into account the flow-sediment interactions;
	A vegetation module was implemented for studying flow-vegetation interactions.

Author(s):
	Gangfeng Ma (gma@odu.edu)
	James T. Kirby (kirby@udel.edu)
	Fengyan Shi (fyshi@udel.edu)

Refernces
	Ma G., Shi F. and Kirby J.T., 2012, Shock-capturing non-hydrostatic model for
	fully dispersive surface wave processes, Ocean Modelling, 43-44, 22-35


