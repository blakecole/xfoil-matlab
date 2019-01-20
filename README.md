# xfoil-matlab
X-FOIL Interface for MATLAB on MacOS

Adapted from Louis Edelman's implementation of Gus Brown's Xfoil Interface.
Modified for MacOS.

System Requirements:
- MacOSX 10.4 or greater with X11 installed.
- XFOIL4MAC (3rd Party Build): http://xfoil4mac.altervista.org/

Improvements:
- Consistent, deterministic output array dimensions: if X-FOIL does not converge, a NaN placeholder is used.
- User-friendly: added additional user warnings and descriptive error messages.
- Updates for modern computers: increase max iterations, number of panels, and panel density near leading and trailing edges.
