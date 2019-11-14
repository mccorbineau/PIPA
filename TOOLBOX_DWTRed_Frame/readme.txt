TOOLBOX_DWTRed_Frame

*****************************************************
* author: Caroline Chaux			    *
* institution: Institut Gaspard Monge - CNRS	    *
*	       Université Paris-Est		    *
* date: Tuesday, Freb., 22nd 2011       	    *
*****************************************************


*****************************************************
* RECOMMENDATIONS:				    *
* This toobox is designed to work with              *
* Matlab 7.0 including the wavelet toolbox 3.0.	    *
*****************************************************

------------------------------------------------------------------------------------
DESCRIPTION:
This toolbox allows to decompose an image via the classical undecimated wavelet transform and a tight frame version.
It consists of 1 subfolder "Wavelab850" that contains a version of the wavelab toolbox.

------------------------------------------------------------------------------------
CAUTION : Add the folder to the path.
	  For that, in Matlab, click on 'File' --> 'Set Path', then 'Add
	  with Subfolders' and choose TOOLBOX_DWTRed_Frame
------------------------------------------------------------------------------------

The filters available are the one available in wavelab.

The files with "adj" suffix designate the adjoint operator.

In the main folder, the m file test_DWT_IT.m is an example script using
the undecimated wavelet transform and demonstrating the perfect reconstruction property.
In the main folder, the m file test_DWT_IT_tight.m is an example script using
the tight frame version of the undecimated wavelet transform and demonstrating the perfect 
reconstruction property.