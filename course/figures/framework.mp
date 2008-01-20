%
% fig2dev (version 3.2.5-alpha5) -L (m)mp version 0.05 --- Preamble
%

%
% mp output driver options:
%  -L mp
%

% +MP-ADDITIONAL-HEADER
verbatimtex
%&latex
\documentclass{article}
\begin{document}
etex
% -MP-ADDITIONAL-HEADER


%SetFigFont macros for latex
verbatimtex
%
\begingroup\makeatletter\ifx\SetFigFont\undefined%
\gdef\SetFigFont#1#2#3#4#5{%
  \reset@font\fontsize{#1}{#2pt}%
  \fontfamily{#3}\fontseries{#4}\fontshape{#5}%
  \selectfont}%
\fi\endgroup%
\ifx\SetFigFontSize\undefined%
\gdef\SetFigFontSize#1#2{%
  \fontsize{#1}{#2pt}%
  \selectfont}%
\fi%
etex

% Now draw the figure
beginfig(0)
% Some reasonable defaults
  labeloffset:=0;
  truecorners:=1;
  bboxmargin:=0;
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 1.89;
  path p;
  p = (170.08, 607.75)
    --(510.24, 451.84);
  draw p withcolor (0.00,0.00,0.00) dashed evenly scaled 1.82;
% End polyline object
% Begin ellipse object
  pickup pencircle scaled 0.94;
  path c;
  c = fullcircle scaled 22.17 yscaled 1.00
         rotated 0.00 shifted (-48.19,349.73);
  draw c withcolor (0.00,0.00,0.00);
% End ellipse object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (-76.54, 267.59)
    --(-48.19, 295.94)
    --(-19.84, 267.59);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (-48.19, 295.94)
    --(-48.19, 298.77)
    --(-48.19, 301.61)
    --(-48.19, 304.44)
    --(-48.19, 307.28)
    --(-48.19, 310.11)
    --(-48.19, 312.94)
    --(-48.19, 315.78)
    --(-48.19, 318.61)
    --(-48.19, 321.45)
    --(-48.19, 324.28)
    --(-48.19, 327.12)
    --(-48.19, 329.95)
    --(-48.19, 332.79)
    --(-48.19, 335.62)
    --(-48.19, 338.46);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (-76.54, 324.28)
    --(-19.84, 324.28);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin ellipse object
  pickup pencircle scaled 0.94;
  path c;
  c = fullcircle scaled 22.17 yscaled 1.00
         rotated 0.00 shifted (-48.19,491.46);
  draw c withcolor (0.00,0.00,0.00);
% End ellipse object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (-76.54, 409.32)
    --(-48.19, 437.67)
    --(-19.84, 409.32);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (-48.19, 437.67)
    --(-48.19, 440.50)
    --(-48.19, 443.34)
    --(-48.19, 446.17)
    --(-48.19, 449.01)
    --(-48.19, 451.84)
    --(-48.19, 454.68)
    --(-48.19, 457.51)
    --(-48.19, 460.35)
    --(-48.19, 463.18)
    --(-48.19, 466.02)
    --(-48.19, 468.85)
    --(-48.19, 471.69)
    --(-48.19, 474.52)
    --(-48.19, 477.35)
    --(-48.19, 480.19);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (-76.54, 466.02)
    --(-19.84, 466.02);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin ellipse object
  pickup pencircle scaled 0.94;
  path c;
  c = fullcircle scaled 22.17 yscaled 1.00
         rotated 0.00 shifted (-48.19,647.37);
  draw c withcolor (0.00,0.00,0.00);
% End ellipse object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (-76.54, 565.23)
    --(-48.19, 593.57)
    --(-19.84, 565.23);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (-48.19, 593.57)
    --(-48.19, 596.41)
    --(-48.19, 599.24)
    --(-48.19, 602.08)
    --(-48.19, 604.91)
    --(-48.19, 607.75)
    --(-48.19, 610.58)
    --(-48.19, 613.42)
    --(-48.19, 616.25)
    --(-48.19, 619.09)
    --(-48.19, 621.92)
    --(-48.19, 624.76)
    --(-48.19, 627.59)
    --(-48.19, 630.43)
    --(-48.19, 633.26)
    --(-48.19, 636.09);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (-76.54, 621.92)
    --(-19.84, 621.92);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt ${\rm S}_1$ etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((598.87,590.55))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt ${\rm S}_2$ etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((597.92,450.71))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt ${\rm S}_3$ etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((597.92,297.64))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin ellipse object
  pickup pencircle scaled 0.94;
  path c;
  c = fullcircle scaled 4.28 yscaled 0.59
         rotated 0.00 shifted (115.78,401.26);
  draw c withcolor (0.00,0.00,0.00);
% End ellipse object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (81.39, 496.94)
    --(80.94, 446.68)
    --(92.85, 446.68)
    --(92.85, 451.46)
    --(143.56, 451.46)
    --(143.56, 446.17)
    --(155.09, 446.17)
    --(156.09, 497.20);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (81.20, 446.43)
    --(75.59, 405.73)
    --(95.12, 402.39)
    --(95.43, 403.21)
    --(145.45, 403.84)
    --(146.96, 403.02)
    --(150.36, 403.09)
    --(164.79, 405.73)
    --(165.61, 406.43)
    --(165.42, 407.37)
    --(154.33, 446.17)
    --(147.53, 446.17)
    --(147.53, 445.35)
    --(87.87, 445.42)
    --(87.87, 446.68)
    --(81.13, 446.68);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (75.72, 405.54)
    --(76.16, 401.95)
    --(95.31, 398.87)
    --(96.13, 399.43)
    --(144.82, 399.87)
    --(146.14, 399.31)
    --(148.54, 399.31)
    --(163.34, 402.58)
    --(165.61, 406.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (155.97, 497.20)
    --(81.57, 496.94)
    --(81.51, 497.07);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (83.84, 494.24)
    --(152.76, 494.49)
    --(152.38, 453.86)
    --(83.34, 453.48)
    --(83.84, 494.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (84.35, 435.46)
    --(84.35, 434.46)
    --(90.08, 434.46)
    --(90.20, 435.40);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (84.60, 435.72)
    --(89.70, 435.72)
    --(89.83, 434.77)
    --(84.47, 434.77)
    --(84.60, 435.78);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (85.35, 443.72)
    --(108.91, 443.65)
    --(108.91, 442.02)
    --(85.29, 442.02)
    --(85.42, 443.78);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (80.13, 432.25)
    --(157.04, 432.19)
    --(163.46, 408.25)
    --(77.35, 407.18)
    --(80.13, 432.19);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (114.27, 403.46)
    --(114.39, 406.24)
    --(128.13, 406.36)
    --(128.25, 403.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (116.16, 405.17)
    --(126.36, 405.17)
    --(126.36, 404.35)
    --(116.16, 404.35)
    --(116.16, 405.17)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (80.94, 433.57)
    --(82.33, 433.57)
    --(82.33, 432.57)
    --(80.94, 432.63)
    --(80.94, 433.57);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (99.59, 411.46)
    --(99.59, 409.20)
    --(124.98, 409.26)
    --(125.10, 411.59)
    --(99.59, 411.46);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (99.59, 409.20)
    --(98.08, 407.56)
    --(98.14, 411.91)
    --(126.43, 412.16)
    --(126.43, 407.75)
    --(125.04, 409.20);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (91.91, 411.46)
    --(91.91, 409.26)
    --(96.94, 409.26)
    --(96.94, 411.53)
    --(91.91, 411.53);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (91.91, 411.46)
    --(90.71, 411.91)
    --(90.71, 407.43)
    --(97.95, 407.43)
    --(97.95, 411.97)
    --(90.71, 411.97);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (90.77, 407.43)
    --(91.91, 409.32);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (97.83, 407.43)
    --(96.94, 409.26);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (97.07, 411.59)
    --(97.89, 411.91);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (85.04, 407.56)
    --(85.04, 411.97)
    --(90.77, 411.97);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (85.04, 407.43)
    --(90.65, 407.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (86.30, 411.46)
    --(85.04, 411.97);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (89.64, 411.46)
    --(90.65, 411.91);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (89.57, 409.45)
    --(90.65, 407.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (86.43, 409.45)
    --(85.10, 407.50);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (85.23, 439.37)
    --(89.89, 439.37)
    --(89.89, 437.48)
    --(85.23, 437.48)
    --(85.23, 439.37)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (86.30, 411.46)
    --(89.64, 411.46)
    --(89.64, 409.45)
    --(86.30, 409.45)
    --(86.30, 411.46)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (79.12, 411.34)
    --(83.97, 411.34)
    --(83.97, 409.07)
    --(79.12, 409.07)
    --(79.12, 411.34)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (78.17, 411.97)
    --(85.04, 411.97)
    --(85.04, 407.50)
    --(78.17, 407.50)
    --(78.17, 411.97)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (78.24, 407.50)
    --(79.12, 409.07);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (79.12, 409.26);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (78.17, 411.97)
    --(79.12, 411.34);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (85.04, 411.97)
    --(83.97, 411.40);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (84.98, 407.50)
    --(83.97, 409.13);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (98.20, 411.78)
    --(99.53, 411.46);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (126.49, 412.09)
    --(131.84, 412.09)
    --(131.84, 407.87)
    --(126.49, 407.87)
    --(126.49, 412.09)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (131.84, 412.09)
    --(138.02, 412.09)
    --(138.02, 407.87)
    --(131.84, 407.87)
    --(131.84, 412.09)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (138.02, 412.16)
    --(142.61, 412.16)
    --(143.43, 408.00);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (143.37, 408.06);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (142.93, 412.16)
    --(143.75, 408.06)
    --(146.39, 408.06)
    --(145.57, 412.09)
    --(142.93, 412.09);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (127.56, 411.46)
    --(131.21, 411.46)
    --(131.21, 409.26)
    --(127.56, 409.26)
    --(127.56, 411.46)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (133.04, 411.53)
    --(137.20, 411.53)
    --(137.20, 409.13)
    --(133.04, 409.13)
    --(133.04, 411.53)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (139.40, 411.28)
    --(142.61, 411.28)
    --(142.61, 408.88)
    --(139.40, 408.88)
    --(139.40, 411.28)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (126.49, 407.75)
    --(127.50, 409.26);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (126.43, 411.97)
    --(127.56, 411.46);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (127.56, 411.53);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (131.84, 412.09)
    --(131.15, 411.53);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (131.84, 407.81)
    --(131.15, 409.26);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (131.91, 407.87)
    --(133.04, 409.13);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (131.84, 412.03)
    --(133.10, 411.53);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (138.02, 412.03)
    --(137.20, 411.53);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (138.02, 407.87)
    --(137.20, 409.13);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (138.14, 407.87)
    --(139.46, 408.94);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (138.02, 412.03)
    --(139.46, 411.28);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (142.61, 412.09)
    --(142.55, 411.21);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (143.24, 408.06)
    --(142.55, 408.94);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (143.50, 416.25)
    --(144.57, 412.35)
    --(150.68, 412.41)
    --(149.54, 416.63)
    --(143.43, 416.50)
    --(143.50, 416.25)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (145.07, 415.69)
    --(145.57, 413.42)
    --(149.67, 413.42)
    --(149.10, 415.87)
    --(145.07, 415.81)
    --(145.07, 415.69)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (145.95, 411.84)
    --(146.83, 408.06)
    --(151.87, 408.13)
    --(150.93, 412.09)
    --(145.95, 411.97)
    --(145.95, 411.84)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (147.78, 411.34)
    --(148.28, 409.13)
    --(151.06, 409.20)
    --(150.61, 411.53)
    --(147.78, 411.46)
    --(147.78, 411.34)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (150.30, 416.63)
    --(151.37, 412.54)
    --(155.97, 412.60)
    --(155.02, 416.69)
    --(150.30, 416.69)
    --(150.30, 416.63)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (151.62, 411.97)
    --(152.50, 408.13)
    --(157.48, 408.25)
    --(156.28, 412.35)
    --(151.56, 412.28)
    --(151.62, 411.97)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (151.87, 415.87)
    --(152.38, 413.73)
    --(155.15, 413.80)
    --(154.71, 416.13)
    --(151.94, 416.13)
    --(151.87, 415.87)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (153.13, 411.53)
    --(153.64, 409.70)
    --(156.54, 409.70)
    --(155.91, 411.84)
    --(153.13, 411.84)
    --(153.13, 411.53)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (156.72, 412.35)
    --(157.86, 408.19)
    --(162.83, 408.25)
    --(161.76, 412.66)
    --(156.72, 412.54)
    --(156.72, 412.35)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (158.24, 411.65)
    --(158.68, 409.76)
    --(161.64, 409.76)
    --(161.26, 411.84)
    --(158.24, 411.78)
    --(158.24, 411.65)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (155.65, 416.57)
    --(156.60, 412.72)
    --(161.57, 412.91)
    --(160.57, 416.76)
    --(155.65, 416.69)
    --(155.65, 416.57)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (156.98, 415.69)
    --(157.42, 414.11)
    --(160.63, 414.11)
    --(160.31, 416.00)
    --(156.98, 415.87)
    --(156.98, 415.69)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (138.08, 416.13)
    --(138.77, 412.54)
    --(144.13, 412.54)
    --(143.24, 416.44)
    --(138.20, 416.50)
    --(138.08, 416.13)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (132.85, 416.50)
    --(133.35, 412.60)
    --(138.39, 412.72)
    --(137.64, 416.50)
    --(132.85, 416.50)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (134.17, 415.69)
    --(134.43, 413.86)
    --(137.32, 413.92)
    --(136.94, 415.87)
    --(134.17, 415.81)
    --(134.17, 415.69)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (139.53, 415.69)
    --(139.84, 414.24)
    --(142.99, 414.30)
    --(142.49, 416.00)
    --(139.40, 416.00)
    --(139.53, 415.69)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (145.45, 420.16)
    --(146.33, 417.07)
    --(154.96, 417.07)
    --(154.02, 420.54)
    --(145.32, 420.72)
    --(145.45, 420.16)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (146.90, 419.97)
    --(147.34, 418.46)
    --(153.95, 418.46)
    --(153.51, 420.03)
    --(146.90, 420.09)
    --(146.90, 419.97)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (154.71, 420.16)
    --(155.46, 417.13)
    --(160.31, 417.32)
    --(159.62, 420.47)
    --(154.71, 420.47)
    --(154.71, 420.16)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (153.70, 424.13)
    --(154.52, 420.66)
    --(159.50, 420.79)
    --(158.49, 424.76)
    --(153.64, 424.63)
    --(153.70, 424.13)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (152.76, 427.78)
    --(153.32, 424.88)
    --(158.43, 424.94)
    --(157.61, 428.03)
    --(152.69, 428.16)
    --(152.76, 427.78)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (148.60, 424.50)
    --(153.13, 424.50)
    --(153.89, 421.04)
    --(149.35, 421.04)
    --(148.60, 424.50)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (143.50, 424.63)
    --(147.91, 424.50)
    --(148.85, 421.29)
    --(144.19, 421.17)
    --(143.43, 424.44)
    --(143.50, 424.63)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (155.78, 419.72)
    --(156.16, 418.20)
    --(159.62, 418.20)
    --(159.24, 419.84)
    --(155.72, 419.78)
    --(155.78, 419.72)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (155.15, 423.69)
    --(155.40, 422.11)
    --(158.74, 422.11)
    --(158.30, 424.06)
    --(155.02, 423.94)
    --(155.15, 423.69)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (154.27, 427.40)
    --(154.46, 426.02)
    --(157.73, 426.02)
    --(157.42, 427.53)
    --(154.27, 427.53)
    --(154.27, 427.40)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (133.48, 412.60)
    --(134.43, 413.86);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (132.85, 416.44)
    --(134.17, 415.81);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (137.64, 416.50)
    --(136.94, 415.81);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (138.46, 412.66)
    --(137.26, 413.92);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (138.90, 412.41)
    --(139.84, 414.36);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (138.14, 416.44)
    --(139.46, 415.94);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (143.31, 416.44)
    --(142.49, 416.00);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (144.06, 412.41)
    --(142.93, 414.30);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (144.69, 412.35)
    --(145.64, 413.42);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (143.43, 416.44)
    --(145.07, 415.81);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (149.48, 416.50)
    --(149.17, 415.81);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (150.74, 412.35)
    --(149.61, 413.48);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (151.50, 412.54)
    --(152.44, 413.80);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (150.43, 416.57)
    --(151.94, 416.00);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (155.15, 416.63)
    --(154.71, 416.13);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (155.91, 412.54)
    --(155.21, 413.86);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (156.66, 412.72)
    --(157.48, 414.17);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (155.72, 416.63)
    --(157.10, 415.87);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (160.63, 416.63)
    --(160.38, 415.94);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (161.64, 412.85)
    --(160.69, 414.17);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (155.46, 417.13)
    --(156.16, 418.27);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (154.77, 420.41)
    --(155.84, 419.78);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (159.69, 420.41)
    --(159.31, 419.78);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (160.38, 417.26)
    --(159.69, 418.27);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (146.83, 408.13)
    --(148.28, 409.13);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (145.89, 411.91)
    --(147.84, 411.34);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (151.94, 408.19)
    --(151.12, 409.13);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (150.99, 411.97)
    --(150.68, 411.53);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (151.56, 412.16)
    --(153.20, 411.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (152.57, 408.13)
    --(153.70, 409.76);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (157.48, 408.19)
    --(156.54, 409.76);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (156.35, 412.22)
    --(155.91, 411.78);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (156.72, 412.47)
    --(158.30, 411.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (157.92, 408.19)
    --(158.74, 409.76);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (162.83, 408.25)
    --(161.64, 409.76);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (161.76, 412.54)
    --(161.32, 411.78);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (144.69, 423.87)
    --(147.65, 423.87)
    --(148.03, 421.98)
    --(145.01, 421.92)
    --(144.63, 423.87)
    --(144.69, 423.87)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (149.67, 423.81)
    --(152.76, 423.87)
    --(153.07, 421.80)
    --(150.24, 421.73)
    --(149.67, 423.69)
    --(149.67, 423.81)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (145.01, 427.84)
    --(145.70, 424.94)
    --(152.88, 424.94)
    --(152.57, 426.83)
    --(152.25, 428.35)
    --(146.20, 428.28)
    --(145.01, 427.97)
    --(145.01, 427.84)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (146.14, 428.22)
    --(146.46, 426.71)
    --(152.50, 426.83);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (139.84, 428.09)
    --(140.47, 424.63)
    --(145.39, 424.76)
    --(144.76, 427.84)
    --(144.25, 428.28)
    --(140.60, 428.22)
    --(139.84, 428.16)
    --(139.84, 428.09)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (140.60, 428.22)
    --(141.04, 426.27)
    --(144.38, 426.46)
    --(144.25, 428.28);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (146.58, 426.65)
    --(145.76, 425.01);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (141.04, 426.27)
    --(140.60, 424.76);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (144.50, 426.33)
    --(145.32, 424.94);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (138.08, 424.31)
    --(138.77, 420.85)
    --(143.56, 420.91)
    --(142.68, 424.44)
    --(138.14, 424.50);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (139.28, 424.06)
    --(139.53, 422.24)
    --(142.55, 422.30)
    --(142.24, 423.94)
    --(139.34, 423.94);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (139.34, 423.81)
    --(138.08, 424.44);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (139.59, 422.24)
    --(138.83, 420.79);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (142.55, 422.24)
    --(143.50, 420.98);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (142.24, 423.94)
    --(142.55, 424.44);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (140.03, 420.41)
    --(140.85, 416.76)
    --(145.89, 416.82)
    --(144.88, 420.54)
    --(140.03, 420.47);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (140.91, 416.76)
    --(141.73, 418.46)
    --(144.63, 418.52)
    --(145.89, 416.88);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (141.73, 418.46)
    --(141.42, 420.09)
    --(144.69, 420.22)
    --(144.69, 418.58);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (141.35, 420.09)
    --(140.22, 420.41);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (144.69, 420.09)
    --(144.88, 420.54);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (134.99, 428.09)
    --(135.56, 424.69)
    --(140.35, 424.69)
    --(139.72, 428.35)
    --(134.93, 428.22);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (135.62, 424.76)
    --(136.50, 425.95)
    --(139.53, 425.95)
    --(140.22, 424.76);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (136.50, 425.95)
    --(136.06, 427.97)
    --(134.99, 428.22);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (136.13, 427.91)
    --(139.34, 427.97)
    --(139.78, 428.22);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (139.34, 427.84)
    --(139.46, 425.95);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (133.10, 424.38)
    --(133.61, 420.60)
    --(138.77, 420.66)
    --(137.83, 424.44)
    --(133.04, 424.50);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (133.67, 420.66)
    --(134.55, 422.36)
    --(137.57, 422.43)
    --(138.71, 420.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (134.61, 422.43)
    --(134.24, 423.94)
    --(133.17, 424.44);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (134.36, 424.00)
    --(137.64, 424.00)
    --(137.89, 424.44);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (137.64, 423.94)
    --(137.64, 422.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (134.87, 420.35)
    --(135.50, 416.76)
    --(140.66, 416.69)
    --(139.72, 420.47)
    --(134.87, 420.47);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (135.50, 416.69)
    --(136.50, 418.27)
    --(139.34, 418.39)
    --(140.72, 416.69);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (136.57, 418.39)
    --(136.50, 419.59)
    --(134.99, 420.41);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (136.57, 419.65)
    --(139.34, 419.65)
    --(139.72, 420.41);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (139.46, 419.59)
    --(139.46, 418.39);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (127.37, 416.13)
    --(127.69, 412.28)
    --(133.17, 412.35)
    --(132.54, 416.19)
    --(127.31, 416.19);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (127.50, 415.69);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (127.69, 412.35)
    --(129.01, 413.92)
    --(131.91, 413.92)
    --(133.17, 412.22);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (129.13, 413.92)
    --(128.82, 415.24)
    --(127.43, 416.19);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (128.82, 415.37)
    --(131.72, 415.43)
    --(132.54, 416.19);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (131.97, 413.92)
    --(131.72, 415.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (131.72, 415.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (129.51, 420.35)
    --(129.83, 416.50)
    --(134.99, 416.76)
    --(134.49, 420.60)
    --(129.51, 420.47);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (129.89, 416.63)
    --(131.21, 417.89)
    --(133.98, 417.89)
    --(134.99, 416.76)
    --(134.87, 416.76);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (129.51, 420.47)
    --(131.15, 419.65)
    --(131.21, 417.89);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (131.21, 419.65)
    --(133.86, 419.78)
    --(134.49, 420.54);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (133.98, 419.72)
    --(133.86, 417.95);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (129.83, 428.35)
    --(130.33, 424.69)
    --(135.12, 424.76)
    --(134.68, 428.35)
    --(129.76, 428.35);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (130.46, 424.63)
    --(131.34, 426.33)
    --(134.11, 426.39)
    --(135.06, 424.69);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (129.83, 428.22)
    --(131.28, 427.46)
    --(131.40, 426.33);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (131.28, 427.46)
    --(134.17, 427.46)
    --(134.17, 426.33);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (134.17, 427.46)
    --(134.74, 428.35);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (127.69, 424.31)
    --(128.19, 420.54)
    --(133.48, 420.66)
    --(132.91, 424.44)
    --(127.69, 424.44);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (128.25, 420.60)
    --(129.32, 422.24)
    --(132.35, 422.24)
    --(133.54, 420.66);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (129.39, 422.24)
    --(129.32, 423.50)
    --(127.69, 424.44);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (129.32, 423.56)
    --(132.22, 423.62)
    --(132.98, 424.44);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (132.22, 423.56)
    --(132.22, 422.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (146.46, 417.01)
    --(147.46, 418.46);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (144.19, 421.10)
    --(145.01, 421.92);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (145.39, 420.60)
    --(146.96, 420.03);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (154.02, 420.60)
    --(153.45, 419.97);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (148.85, 421.23)
    --(148.03, 421.98);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (147.97, 424.57)
    --(147.72, 423.81);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (143.50, 424.50)
    --(144.69, 423.81);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (148.66, 424.57)
    --(149.67, 423.81);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (149.35, 421.04)
    --(150.24, 421.73);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (154.02, 421.04)
    --(153.07, 421.73);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (153.13, 424.57)
    --(152.82, 423.81);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (152.94, 423.81);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (153.39, 424.88)
    --(154.46, 426.02);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (152.82, 428.09)
    --(154.39, 427.46);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (157.67, 428.03)
    --(157.42, 427.53);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (158.49, 424.88)
    --(157.73, 426.02);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (154.90, 417.20)
    --(153.95, 418.46);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (154.71, 420.66)
    --(155.53, 422.11);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (153.76, 424.57)
    --(155.21, 423.87);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (158.61, 424.69)
    --(158.36, 423.94);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (159.56, 420.79)
    --(158.74, 422.11);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (124.28, 420.35)
    --(124.54, 416.63)
    --(129.76, 416.50)
    --(129.26, 420.47)
    --(124.28, 420.47);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (124.66, 416.63)
    --(125.54, 417.89)
    --(128.69, 417.95)
    --(129.70, 416.57);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (125.54, 417.89)
    --(125.54, 419.65)
    --(124.47, 420.47);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (125.54, 419.59)
    --(128.63, 419.59)
    --(129.20, 420.47);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (128.63, 419.59)
    --(128.69, 417.89);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (121.76, 416.31)
    --(122.14, 412.09)
    --(127.50, 412.22)
    --(127.18, 416.25)
    --(121.76, 416.31);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (122.08, 412.16)
    --(123.15, 413.48)
    --(126.36, 413.42)
    --(127.50, 412.22);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (126.43, 413.48)
    --(126.36, 415.37)
    --(127.18, 416.19);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (126.36, 415.43)
    --(122.96, 415.37)
    --(123.09, 413.48);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (123.02, 415.31)
    --(121.83, 416.31);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (118.93, 420.47)
    --(119.12, 416.57)
    --(124.41, 416.57)
    --(123.97, 420.47)
    --(118.93, 420.47);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (119.12, 416.63)
    --(120.19, 417.89)
    --(123.09, 417.95)
    --(124.47, 416.63);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (120.25, 417.89)
    --(120.19, 419.65)
    --(118.93, 420.47);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (120.31, 419.65)
    --(123.09, 419.72)
    --(123.97, 420.41);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (123.09, 419.65)
    --(123.09, 417.89);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (116.41, 416.31)
    --(116.41, 412.28)
    --(121.89, 412.28)
    --(121.57, 416.25)
    --(116.41, 416.31);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (116.41, 412.28)
    --(117.61, 413.67)
    --(120.82, 413.67)
    --(121.76, 412.28);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (117.61, 413.73)
    --(117.48, 415.24)
    --(116.41, 416.31);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (117.42, 415.37)
    --(120.88, 415.31)
    --(121.57, 416.19);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (120.88, 415.31)
    --(120.82, 413.67);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (122.46, 424.38)
    --(122.90, 420.66)
    --(127.94, 420.66)
    --(127.37, 424.44)
    --(122.46, 424.38);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (123.02, 420.66)
    --(123.84, 422.05)
    --(126.93, 422.05)
    --(127.94, 420.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (123.84, 422.11)
    --(123.72, 423.62)
    --(122.58, 424.38);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (123.78, 423.69)
    --(127.06, 423.75)
    --(127.37, 424.44);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (127.06, 423.62)
    --(126.93, 422.11);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (124.66, 428.16)
    --(125.04, 424.76)
    --(129.89, 424.76)
    --(129.39, 428.35)
    --(124.66, 428.35)
    --(124.66, 428.16);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (125.17, 424.76)
    --(125.92, 425.95)
    --(129.01, 425.95)
    --(129.95, 424.82);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (129.07, 425.95)
    --(128.88, 427.59)
    --(129.32, 428.35);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (128.82, 427.59)
    --(125.61, 427.59)
    --(125.92, 425.95);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (125.54, 427.65)
    --(124.66, 428.35);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (119.50, 428.60)
    --(119.81, 424.82)
    --(124.66, 424.82)
    --(124.09, 428.54)
    --(119.50, 428.54);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (119.81, 424.82)
    --(120.82, 425.95)
    --(123.59, 426.02)
    --(124.66, 424.88);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (119.50, 428.47)
    --(120.57, 427.97)
    --(120.82, 425.89);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (120.69, 427.97)
    --(123.59, 427.97)
    --(124.16, 428.60);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (123.72, 427.97)
    --(123.65, 425.95)
    --(123.65, 426.02);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (117.10, 424.44)
    --(117.35, 420.66)
    --(122.71, 420.66)
    --(122.08, 424.44)
    --(117.17, 424.38);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (117.35, 420.66)
    --(118.43, 421.98)
    --(121.57, 421.98)
    --(122.71, 420.79);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (118.49, 421.98)
    --(118.36, 423.43)
    --(117.23, 424.38);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (118.43, 423.43)
    --(121.57, 423.43)
    --(122.08, 424.44);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (121.57, 423.37)
    --(121.51, 421.98);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (114.39, 428.22)
    --(114.46, 424.82)
    --(119.50, 424.94)
    --(119.12, 428.41)
    --(114.39, 428.35);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (114.46, 424.88)
    --(115.28, 425.89)
    --(118.68, 425.95)
    --(119.43, 424.94);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (115.34, 425.83)
    --(115.21, 427.78)
    --(114.39, 428.41);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (115.34, 427.78)
    --(118.68, 427.84)
    --(119.12, 428.35);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (118.68, 427.72)
    --(118.74, 425.95);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (111.94, 424.38)
    --(112.06, 420.72)
    --(117.04, 420.72)
    --(116.72, 424.44)
    --(111.94, 424.50);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (112.06, 420.72)
    --(112.94, 421.67)
    --(116.22, 421.73)
    --(117.04, 420.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (112.94, 421.73)
    --(112.76, 423.62)
    --(112.06, 424.44);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (112.76, 423.62)
    --(116.22, 423.62)
    --(116.72, 424.44);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (116.22, 423.62)
    --(116.22, 421.73);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (113.26, 420.41)
    --(113.20, 416.88)
    --(118.80, 416.76)
    --(118.49, 420.47)
    --(113.26, 420.54);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (113.26, 416.76)
    --(114.27, 417.70)
    --(117.92, 417.70)
    --(118.74, 416.82);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (114.27, 417.70)
    --(114.27, 419.65)
    --(113.26, 420.41);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (114.08, 419.78)
    --(117.98, 419.78)
    --(118.43, 420.41);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (118.49, 420.35);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (118.05, 419.65)
    --(117.98, 417.70);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (110.55, 416.31)
    --(110.80, 412.28)
    --(116.03, 412.35)
    --(115.72, 416.44)
    --(110.55, 416.44);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (110.87, 412.35)
    --(111.81, 413.42)
    --(115.21, 413.42)
    --(115.91, 412.41);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (111.87, 413.42)
    --(111.75, 415.50)
    --(110.55, 416.44);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (111.81, 415.50)
    --(115.09, 415.50)
    --(115.72, 416.38);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (115.09, 415.50)
    --(115.15, 413.42);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (109.10, 428.22)
    --(109.35, 425.01)
    --(114.14, 425.01)
    --(113.70, 428.35)
    --(109.04, 428.35);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (109.35, 425.13)
    --(110.17, 425.95)
    --(113.32, 425.95)
    --(114.08, 425.07);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (110.17, 425.95)
    --(110.11, 427.53)
    --(109.10, 428.28);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (110.17, 427.59)
    --(113.26, 427.59)
    --(113.64, 428.28);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (113.32, 427.53)
    --(113.32, 425.95);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (106.52, 424.50)
    --(106.65, 420.85)
    --(111.81, 420.72)
    --(111.43, 424.44)
    --(106.46, 424.50);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (106.71, 420.85)
    --(107.59, 422.11)
    --(110.74, 422.17)
    --(111.81, 420.79);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (107.72, 422.24)
    --(107.65, 423.62)
    --(106.52, 424.50);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (107.72, 423.62)
    --(110.80, 423.62)
    --(111.43, 424.44);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (110.80, 423.56)
    --(110.80, 422.17);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (107.78, 420.41)
    --(107.91, 416.69)
    --(112.94, 416.76)
    --(112.82, 420.35)
    --(107.78, 420.47);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (107.91, 416.63)
    --(108.85, 417.64)
    --(112.25, 417.64)
    --(112.94, 416.76);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (108.91, 417.70)
    --(108.91, 419.65)
    --(107.78, 420.47);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (108.98, 419.65)
    --(112.38, 419.65)
    --(112.88, 420.35);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (112.38, 419.65)
    --(112.25, 417.64);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (105.13, 416.31)
    --(105.26, 412.22)
    --(110.49, 412.28)
    --(110.30, 416.31)
    --(105.13, 416.31);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (105.26, 412.22)
    --(106.39, 413.48)
    --(109.61, 413.54)
    --(110.49, 412.28);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (106.52, 413.48)
    --(106.39, 415.43)
    --(105.13, 416.31);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (106.39, 415.50)
    --(109.61, 415.50)
    --(110.30, 416.38);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (109.61, 415.43)
    --(109.61, 413.48);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (104.00, 428.47)
    --(104.00, 424.76)
    --(108.72, 424.76)
    --(108.41, 428.54)
    --(104.06, 428.54);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (103.94, 424.82)
    --(104.88, 425.89)
    --(107.78, 425.89)
    --(108.72, 424.82);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (105.01, 425.89)
    --(105.01, 427.84)
    --(104.00, 428.60);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (105.01, 427.91)
    --(107.84, 427.91)
    --(108.47, 428.54);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (107.84, 427.91)
    --(107.78, 425.83);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (101.35, 424.50)
    --(101.29, 420.79)
    --(106.27, 420.79)
    --(106.14, 424.50)
    --(101.35, 424.57);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (101.29, 420.79)
    --(102.30, 421.98)
    --(105.39, 421.98)
    --(106.27, 420.85);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (102.36, 421.98)
    --(102.43, 423.75)
    --(101.35, 424.50);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (102.55, 423.87)
    --(105.45, 423.87)
    --(106.14, 424.50);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (105.51, 423.75)
    --(105.45, 421.98);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (102.43, 420.54)
    --(102.43, 416.63)
    --(107.40, 416.63)
    --(107.40, 420.47)
    --(102.43, 420.54);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (102.43, 416.69)
    --(103.31, 417.70)
    --(106.71, 417.70)
    --(107.40, 416.69);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (103.43, 417.70)
    --(103.37, 419.53)
    --(102.43, 420.47);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (103.43, 419.53)
    --(106.71, 419.53)
    --(107.40, 420.47);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (106.71, 419.53)
    --(106.71, 417.70);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (99.65, 416.31)
    --(99.65, 412.22)
    --(104.82, 412.22)
    --(104.82, 416.31)
    --(99.65, 416.31)
    --(99.72, 416.38);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (99.72, 412.35)
    --(100.72, 413.92)
    --(104.00, 413.92)
    --(104.82, 412.22);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (100.79, 413.92)
    --(100.85, 415.50)
    --(99.72, 416.31);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (100.85, 415.50)
    --(103.81, 415.50)
    --(104.76, 416.19);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (103.87, 415.50)
    --(103.94, 413.92);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (98.27, 424.82)
    --(103.62, 424.82)
    --(103.43, 428.41)
    --(98.39, 428.47)
    --(98.46, 424.94);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (98.46, 424.82)
    --(99.72, 426.20)
    --(102.55, 426.27)
    --(103.62, 424.94);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (99.72, 426.20)
    --(99.65, 428.22)
    --(98.39, 428.47);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (99.72, 428.16)
    --(102.87, 428.09)
    --(103.50, 428.41);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (103.50, 428.22);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (102.93, 428.03)
    --(102.55, 426.20);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (93.92, 416.44)
    --(94.05, 412.22)
    --(99.21, 412.22)
    --(99.21, 416.44)
    --(93.92, 416.57);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (94.11, 412.22)
    --(95.12, 413.48)
    --(98.08, 413.48)
    --(99.28, 412.28);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (98.08, 413.48)
    --(98.08, 415.87)
    --(99.28, 416.44);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (97.95, 415.94)
    --(94.99, 415.87)
    --(94.05, 416.44);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (95.06, 415.81)
    --(95.24, 413.42);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (96.88, 420.35)
    --(96.94, 416.76)
    --(102.17, 416.76)
    --(101.86, 420.47)
    --(96.88, 420.54);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (97.01, 416.76)
    --(97.83, 417.89)
    --(101.29, 417.89)
    --(102.17, 416.76);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (101.29, 417.95)
    --(101.29, 420.22)
    --(101.86, 420.47);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (97.83, 417.95)
    --(97.76, 420.22)
    --(96.94, 420.47);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (97.89, 420.16)
    --(101.29, 420.16);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (95.94, 424.25)
    --(95.94, 420.98)
    --(100.98, 420.98)
    --(100.98, 424.44)
    --(95.87, 424.44);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (95.87, 420.98)
    --(96.63, 422.17)
    --(100.16, 422.17)
    --(101.04, 420.98);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (100.16, 422.24)
    --(100.16, 424.06)
    --(100.98, 424.44);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (100.16, 424.06)
    --(96.69, 424.13)
    --(95.81, 424.44);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (96.63, 422.24)
    --(96.63, 424.13);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (93.29, 428.60)
    --(93.29, 424.82)
    --(98.27, 424.76)
    --(98.08, 428.54)
    --(93.42, 428.54);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (93.29, 424.82)
    --(94.30, 426.20)
    --(97.20, 426.20)
    --(98.20, 424.82);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (97.20, 426.20)
    --(97.26, 428.28)
    --(98.14, 428.54);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (97.26, 428.28)
    --(94.36, 428.28)
    --(93.29, 428.54);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (94.24, 428.28)
    --(94.36, 426.20);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (88.31, 428.47)
    --(88.19, 424.88)
    --(92.91, 424.88)
    --(92.85, 428.54)
    --(88.31, 428.54);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (88.31, 424.94)
    --(89.07, 426.20)
    --(92.09, 426.27)
    --(92.91, 424.88);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (92.09, 426.27)
    --(92.09, 428.28)
    --(92.91, 428.54);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (92.09, 428.28)
    --(89.20, 428.28)
    --(89.20, 426.27);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (89.20, 428.28)
    --(88.25, 428.54);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (90.52, 424.44)
    --(90.39, 420.60)
    --(95.43, 420.60)
    --(95.56, 424.50)
    --(90.52, 424.63);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (90.33, 420.66)
    --(91.40, 422.05)
    --(94.55, 422.05)
    --(95.37, 420.60);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (94.61, 422.05)
    --(94.61, 424.00)
    --(95.56, 424.50);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (94.55, 424.00)
    --(91.28, 424.00)
    --(91.34, 422.11);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (91.34, 423.94)
    --(90.52, 424.57);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (91.59, 420.35)
    --(91.59, 416.76)
    --(96.63, 416.76)
    --(96.63, 420.47)
    --(91.59, 420.47);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (91.65, 416.82)
    --(92.66, 418.02)
    --(95.87, 418.02)
    --(96.57, 416.88);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (95.94, 418.02)
    --(95.87, 420.16)
    --(96.57, 420.47);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (96.00, 420.16)
    --(92.41, 420.22)
    --(92.66, 417.95);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (92.47, 420.22)
    --(91.59, 420.47);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (88.57, 416.06)
    --(88.57, 412.41)
    --(93.67, 412.41)
    --(93.67, 416.25)
    --(88.57, 416.38)
    --(88.57, 416.00);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (88.57, 415.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (88.57, 412.41)
    --(89.51, 413.73)
    --(92.79, 413.73)
    --(93.61, 412.47);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (93.42, 412.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (92.79, 413.67)
    --(92.91, 415.87)
    --(93.73, 416.31);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (92.91, 415.87)
    --(89.32, 415.87)
    --(89.45, 413.80);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (89.45, 415.81)
    --(88.57, 416.38);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (85.92, 420.28)
    --(85.98, 416.63)
    --(91.40, 416.57)
    --(91.15, 420.41)
    --(85.92, 420.41);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (86.05, 416.57)
    --(86.93, 418.02)
    --(90.39, 418.02)
    --(91.40, 416.63);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (90.46, 418.02)
    --(90.39, 419.97)
    --(91.21, 420.41);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (90.46, 419.97)
    --(86.74, 419.91)
    --(86.93, 418.08);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (86.80, 419.91)
    --(85.92, 420.35);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (85.10, 424.25)
    --(85.04, 420.54)
    --(90.33, 420.60)
    --(90.27, 424.57)
    --(85.10, 424.57)
    --(85.10, 424.25);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (85.17, 420.60)
    --(85.92, 421.98)
    --(89.13, 421.98)
    --(90.39, 420.60);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (89.20, 421.98)
    --(89.13, 424.00)
    --(90.33, 424.50);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (89.20, 423.94)
    --(85.92, 423.94)
    --(85.98, 421.98);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (86.05, 424.00)
    --(85.10, 424.63);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (82.96, 428.47)
    --(82.77, 424.94)
    --(87.81, 424.88)
    --(87.69, 428.47)
    --(82.96, 428.54);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (82.83, 425.01)
    --(83.84, 426.20)
    --(86.61, 426.20)
    --(87.81, 424.94);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (86.55, 426.20)
    --(86.61, 427.97)
    --(87.75, 428.54);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (86.61, 427.97)
    --(83.97, 427.97)
    --(83.91, 426.20);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (83.97, 427.91)
    --(83.02, 428.47);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (83.02, 428.60);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (78.87, 416.06)
    --(78.30, 412.22)
    --(88.06, 412.41)
    --(88.25, 416.19)
    --(78.80, 416.19);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (78.36, 412.28)
    --(79.37, 413.67)
    --(86.43, 413.80)
    --(88.13, 412.41);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (86.49, 413.73)
    --(86.74, 415.56)
    --(88.19, 416.19);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (88.25, 416.19);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (86.74, 415.56)
    --(79.50, 415.62)
    --(79.37, 413.67);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (79.50, 415.56)
    --(78.93, 416.13);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (79.18, 420.35)
    --(78.93, 416.50)
    --(85.73, 416.50)
    --(85.54, 420.35)
    --(79.18, 420.35);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (79.12, 416.57)
    --(80.19, 417.64)
    --(84.72, 417.57)
    --(85.73, 416.63);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (84.79, 417.51)
    --(84.79, 419.65)
    --(85.54, 420.35);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (84.72, 419.59)
    --(80.13, 419.59)
    --(80.19, 417.51);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (80.19, 419.59)
    --(79.31, 420.35);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (79.43, 420.72)
    --(84.91, 420.72)
    --(84.72, 424.50)
    --(80.13, 424.57)
    --(79.31, 420.66);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (79.50, 420.72)
    --(80.50, 422.43)
    --(82.90, 422.49)
    --(84.98, 420.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (80.44, 422.43)
    --(80.82, 424.06)
    --(80.13, 424.57);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (80.82, 424.06)
    --(82.90, 424.06)
    --(82.90, 422.49);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (82.90, 424.13)
    --(84.60, 424.44);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (80.31, 428.85)
    --(82.83, 428.85)
    --(82.58, 424.69)
    --(80.06, 424.69)
    --(80.38, 428.85);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (80.82, 431.56)
    --(80.50, 429.10)
    --(85.48, 428.98)
    --(85.86, 431.43)
    --(80.82, 431.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (80.57, 429.23)
    --(81.13, 430.36)
    --(84.72, 430.36)
    --(85.48, 429.04);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (84.79, 430.36)
    --(84.85, 431.37);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (81.20, 430.43)
    --(81.51, 431.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (86.68, 431.37)
    --(86.49, 428.98)
    --(91.65, 428.98)
    --(91.97, 431.24)
    --(86.80, 431.31);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (86.43, 428.91)
    --(87.24, 429.67)
    --(91.09, 429.80)
    --(91.78, 428.98);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (87.24, 429.80)
    --(87.24, 431.18)
    --(87.31, 431.18);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (91.09, 429.67)
    --(91.21, 431.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (92.47, 431.12)
    --(92.41, 429.04)
    --(96.44, 428.98)
    --(96.63, 431.24)
    --(92.41, 431.18);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (92.41, 428.98)
    --(93.23, 429.86)
    --(95.56, 429.86)
    --(96.44, 428.98);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (95.62, 429.92)
    --(95.69, 430.99);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (93.23, 429.92)
    --(93.29, 431.18);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (97.07, 431.12)
    --(97.07, 428.79)
    --(101.29, 428.79)
    --(101.35, 431.12)
    --(97.13, 431.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (97.20, 428.85)
    --(97.70, 429.67)
    --(100.54, 429.67)
    --(101.42, 428.79);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (97.76, 429.67)
    --(97.83, 431.18);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (100.54, 429.67)
    --(100.54, 431.06);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (101.86, 431.12)
    --(101.80, 428.79)
    --(106.08, 428.79)
    --(106.08, 431.24)
    --(101.86, 431.31);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (101.98, 431.31);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (101.86, 428.79)
    --(102.74, 429.80)
    --(105.39, 429.73)
    --(106.08, 428.85);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (102.74, 429.73)
    --(102.68, 431.31);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (105.32, 429.73)
    --(105.39, 431.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (107.02, 431.18)
    --(106.96, 428.79)
    --(111.06, 428.79)
    --(111.06, 431.24)
    --(106.96, 431.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (106.96, 428.79)
    --(107.53, 429.80)
    --(110.24, 429.80)
    --(111.12, 428.85);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (107.53, 429.80)
    --(107.53, 431.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (110.17, 429.80)
    --(110.17, 431.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (111.37, 431.18)
    --(111.37, 428.72)
    --(115.59, 428.72)
    --(115.59, 431.18)
    --(111.43, 431.18);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (111.37, 431.31);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (111.37, 428.66)
    --(111.87, 429.48)
    --(114.90, 429.48)
    --(115.59, 428.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (114.90, 429.54)
    --(114.90, 431.18);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (111.87, 429.48)
    --(111.87, 431.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (115.97, 431.12)
    --(115.97, 428.66)
    --(119.94, 428.72)
    --(119.87, 431.24)
    --(115.91, 431.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (115.97, 428.66)
    --(116.54, 429.48)
    --(119.18, 429.48)
    --(119.94, 428.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (116.60, 429.48)
    --(116.60, 431.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (119.18, 429.48)
    --(119.18, 431.18);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (120.38, 431.18)
    --(120.31, 428.79)
    --(124.60, 428.72)
    --(124.60, 431.18)
    --(120.38, 431.18);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (120.44, 428.72)
    --(121.07, 429.48)
    --(124.09, 429.48)
    --(124.54, 428.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (125.10, 431.18)
    --(125.10, 428.66)
    --(129.51, 428.66)
    --(129.51, 431.24)
    --(125.10, 431.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (125.10, 428.66)
    --(125.61, 429.48)
    --(129.13, 429.48)
    --(129.51, 428.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (125.67, 429.48)
    --(125.61, 431.31);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (129.13, 429.48)
    --(129.13, 431.31);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (121.07, 429.48)
    --(121.07, 431.18);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (124.22, 429.48)
    --(124.22, 431.18);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (129.95, 431.24)
    --(129.95, 428.66)
    --(134.43, 428.72)
    --(134.43, 431.24)
    --(129.95, 431.37);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (129.95, 428.72)
    --(130.65, 429.48)
    --(133.73, 429.48)
    --(134.43, 428.79);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (130.77, 429.48)
    --(130.77, 431.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (133.80, 429.48)
    --(133.80, 431.31);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (134.74, 431.18)
    --(134.74, 428.66)
    --(139.40, 428.66)
    --(139.28, 431.31)
    --(134.68, 431.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (134.80, 428.66)
    --(135.37, 429.35)
    --(138.65, 429.35)
    --(139.46, 428.66);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (135.50, 429.42)
    --(135.50, 431.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (138.65, 429.35)
    --(138.65, 431.31);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (139.65, 431.37)
    --(139.65, 428.66)
    --(143.87, 428.66)
    --(143.75, 431.37)
    --(139.65, 431.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (139.72, 428.66)
    --(140.47, 429.48)
    --(143.18, 429.48)
    --(143.87, 428.66);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (143.18, 429.48)
    --(143.18, 431.37);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (140.54, 429.48)
    --(140.47, 431.37);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (144.31, 431.37)
    --(144.31, 428.66)
    --(148.54, 428.66)
    --(148.54, 431.37)
    --(144.31, 431.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (144.31, 428.60)
    --(145.13, 429.54)
    --(147.91, 429.54)
    --(148.54, 428.66);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (145.20, 429.48)
    --(145.20, 431.37);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (147.91, 429.48)
    --(147.91, 431.37);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (148.85, 431.31)
    --(148.79, 428.60)
    --(152.57, 428.60)
    --(152.57, 431.24)
    --(148.85, 431.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (148.85, 428.66)
    --(149.54, 429.54)
    --(152.13, 429.54)
    --(152.50, 428.66);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (149.67, 429.54)
    --(149.67, 431.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (152.25, 429.48)
    --(152.25, 431.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (153.01, 431.12)
    --(153.01, 428.54)
    --(157.61, 428.54)
    --(156.98, 431.31)
    --(153.01, 431.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (153.13, 428.54)
    --(153.83, 429.42)
    --(157.04, 429.42)
    --(157.61, 428.60);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (157.10, 429.42)
    --(156.60, 431.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (153.95, 429.48)
    --(153.76, 431.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (155.53, 433.57)
    --(157.17, 433.57)
    --(157.17, 432.76)
    --(155.53, 432.76)
    --(155.53, 433.57)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (107.21, 497.01)
    --(107.46, 495.06)
    --(130.52, 495.24)
    --(130.58, 497.13);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (125.04, 411.53)
    --(126.49, 412.09);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin ellipse object
  pickup pencircle scaled 0.94;
  path c;
  c = fullcircle scaled 4.28 yscaled 0.59
         rotated 0.00 shifted (115.78,552.50);
  draw c withcolor (0.00,0.00,0.00);
% End ellipse object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (81.39, 648.19)
    --(80.94, 597.92)
    --(92.85, 597.92)
    --(92.85, 602.71)
    --(143.56, 602.71)
    --(143.56, 597.42)
    --(155.09, 597.42)
    --(156.09, 648.44);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (81.20, 597.67)
    --(75.59, 556.98)
    --(95.12, 553.64)
    --(95.43, 554.46)
    --(145.45, 555.09)
    --(146.96, 554.27)
    --(150.36, 554.33)
    --(164.79, 556.98)
    --(165.61, 557.67)
    --(165.42, 558.61)
    --(154.33, 597.42)
    --(147.53, 597.42)
    --(147.53, 596.60)
    --(87.87, 596.66)
    --(87.87, 597.92)
    --(81.13, 597.92);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (75.72, 556.79)
    --(76.16, 553.20)
    --(95.31, 550.11)
    --(96.13, 550.68)
    --(144.82, 551.12)
    --(146.14, 550.55)
    --(148.54, 550.55)
    --(163.34, 553.83)
    --(165.61, 557.48);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (155.97, 648.44)
    --(81.57, 648.19)
    --(81.51, 648.31);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (83.84, 645.48)
    --(152.76, 645.73)
    --(152.38, 605.10)
    --(83.34, 604.72)
    --(83.84, 645.67);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (84.35, 586.71)
    --(84.35, 585.70)
    --(90.08, 585.70)
    --(90.20, 586.65);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (84.60, 586.96)
    --(89.70, 586.96)
    --(89.83, 586.02)
    --(84.47, 586.02)
    --(84.60, 587.02);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (85.35, 594.96)
    --(108.91, 594.90)
    --(108.91, 593.26)
    --(85.29, 593.26)
    --(85.42, 595.02);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (80.13, 583.50)
    --(157.04, 583.43)
    --(163.46, 559.50)
    --(77.35, 558.43)
    --(80.13, 583.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (114.27, 554.71)
    --(114.39, 557.48)
    --(128.13, 557.61)
    --(128.25, 554.96);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (116.16, 556.41)
    --(126.36, 556.41)
    --(126.36, 555.59)
    --(116.16, 555.59)
    --(116.16, 556.41)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (80.94, 584.82)
    --(82.33, 584.82)
    --(82.33, 583.81)
    --(80.94, 583.87)
    --(80.94, 584.82);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (99.59, 562.71)
    --(99.59, 560.44)
    --(124.98, 560.50)
    --(125.10, 562.83)
    --(99.59, 562.71);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (99.59, 560.44)
    --(98.08, 558.80)
    --(98.14, 563.15)
    --(126.43, 563.40)
    --(126.43, 558.99)
    --(125.04, 560.44);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (91.91, 562.71)
    --(91.91, 560.50)
    --(96.94, 560.50)
    --(96.94, 562.77)
    --(91.91, 562.77);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (91.91, 562.71)
    --(90.71, 563.15)
    --(90.71, 558.68)
    --(97.95, 558.68)
    --(97.95, 563.21)
    --(90.71, 563.21);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (90.77, 558.68)
    --(91.91, 560.57);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (97.83, 558.68)
    --(96.94, 560.50);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (97.07, 562.83)
    --(97.89, 563.15);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (85.04, 558.80)
    --(85.04, 563.21)
    --(90.77, 563.21);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (85.04, 558.68)
    --(90.65, 558.68);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (86.30, 562.71)
    --(85.04, 563.21);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (89.64, 562.71)
    --(90.65, 563.15);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (89.57, 560.69)
    --(90.65, 558.68);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (86.43, 560.69)
    --(85.10, 558.74);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (85.23, 590.61)
    --(89.89, 590.61)
    --(89.89, 588.72)
    --(85.23, 588.72)
    --(85.23, 590.61)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (86.30, 562.71)
    --(89.64, 562.71)
    --(89.64, 560.69)
    --(86.30, 560.69)
    --(86.30, 562.71)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (79.12, 562.58)
    --(83.97, 562.58)
    --(83.97, 560.31)
    --(79.12, 560.31)
    --(79.12, 562.58)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (78.17, 563.21)
    --(85.04, 563.21)
    --(85.04, 558.74)
    --(78.17, 558.74)
    --(78.17, 563.21)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (78.24, 558.74)
    --(79.12, 560.31);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (79.12, 560.50);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (78.17, 563.21)
    --(79.12, 562.58);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (85.04, 563.21)
    --(83.97, 562.65);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (84.98, 558.74)
    --(83.97, 560.38);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (98.20, 563.02)
    --(99.53, 562.71);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (126.49, 563.34)
    --(131.84, 563.34)
    --(131.84, 559.12)
    --(126.49, 559.12)
    --(126.49, 563.34)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (131.84, 563.34)
    --(138.02, 563.34)
    --(138.02, 559.12)
    --(131.84, 559.12)
    --(131.84, 563.34)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (138.02, 563.40)
    --(142.61, 563.40)
    --(143.43, 559.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (143.37, 559.31);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (142.93, 563.40)
    --(143.75, 559.31)
    --(146.39, 559.31)
    --(145.57, 563.34)
    --(142.93, 563.34);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (127.56, 562.71)
    --(131.21, 562.71)
    --(131.21, 560.50)
    --(127.56, 560.50)
    --(127.56, 562.71)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (133.04, 562.77)
    --(137.20, 562.77)
    --(137.20, 560.38)
    --(133.04, 560.38)
    --(133.04, 562.77)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (139.40, 562.52)
    --(142.61, 562.52)
    --(142.61, 560.13)
    --(139.40, 560.13)
    --(139.40, 562.52)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (126.49, 558.99)
    --(127.50, 560.50);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (126.43, 563.21)
    --(127.56, 562.71);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (127.56, 562.77);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (131.84, 563.34)
    --(131.15, 562.77);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (131.84, 559.06)
    --(131.15, 560.50);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (131.91, 559.12)
    --(133.04, 560.38);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (131.84, 563.28)
    --(133.10, 562.77);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (138.02, 563.28)
    --(137.20, 562.77);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (138.02, 559.12)
    --(137.20, 560.38);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (138.14, 559.12)
    --(139.46, 560.19);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (138.02, 563.28)
    --(139.46, 562.52);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (142.61, 563.34)
    --(142.55, 562.46);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (143.24, 559.31)
    --(142.55, 560.19);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (143.50, 567.50)
    --(144.57, 563.59)
    --(150.68, 563.65)
    --(149.54, 567.87)
    --(143.43, 567.75)
    --(143.50, 567.50)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (145.07, 566.93)
    --(145.57, 564.66)
    --(149.67, 564.66)
    --(149.10, 567.12)
    --(145.07, 567.06)
    --(145.07, 566.93)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (145.95, 563.09)
    --(146.83, 559.31)
    --(151.87, 559.37)
    --(150.93, 563.34)
    --(145.95, 563.21)
    --(145.95, 563.09)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (147.78, 562.58)
    --(148.28, 560.38)
    --(151.06, 560.44)
    --(150.61, 562.77)
    --(147.78, 562.71)
    --(147.78, 562.58)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (150.30, 567.87)
    --(151.37, 563.78)
    --(155.97, 563.84)
    --(155.02, 567.94)
    --(150.30, 567.94)
    --(150.30, 567.87)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (151.62, 563.21)
    --(152.50, 559.37)
    --(157.48, 559.50)
    --(156.28, 563.59)
    --(151.56, 563.53)
    --(151.62, 563.21)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (151.87, 567.12)
    --(152.38, 564.98)
    --(155.15, 565.04)
    --(154.71, 567.37)
    --(151.94, 567.37)
    --(151.87, 567.12)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (153.13, 562.77)
    --(153.64, 560.94)
    --(156.54, 560.94)
    --(155.91, 563.09)
    --(153.13, 563.09)
    --(153.13, 562.77)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (156.72, 563.59)
    --(157.86, 559.43)
    --(162.83, 559.50)
    --(161.76, 563.91)
    --(156.72, 563.78)
    --(156.72, 563.59)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (158.24, 562.90)
    --(158.68, 561.01)
    --(161.64, 561.01)
    --(161.26, 563.09)
    --(158.24, 563.02)
    --(158.24, 562.90)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (155.65, 567.81)
    --(156.60, 563.97)
    --(161.57, 564.16)
    --(160.57, 568.00)
    --(155.65, 567.94)
    --(155.65, 567.81)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (156.98, 566.93)
    --(157.42, 565.35)
    --(160.63, 565.35)
    --(160.31, 567.24)
    --(156.98, 567.12)
    --(156.98, 566.93)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (138.08, 567.37)
    --(138.77, 563.78)
    --(144.13, 563.78)
    --(143.24, 567.69)
    --(138.20, 567.75)
    --(138.08, 567.37)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (132.85, 567.75)
    --(133.35, 563.84)
    --(138.39, 563.97)
    --(137.64, 567.75)
    --(132.85, 567.75)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (134.17, 566.93)
    --(134.43, 565.10)
    --(137.32, 565.17)
    --(136.94, 567.12)
    --(134.17, 567.06)
    --(134.17, 566.93)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (139.53, 566.93)
    --(139.84, 565.48)
    --(142.99, 565.54)
    --(142.49, 567.24)
    --(139.40, 567.24)
    --(139.53, 566.93)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (145.45, 571.40)
    --(146.33, 568.31)
    --(154.96, 568.31)
    --(154.02, 571.78)
    --(145.32, 571.97)
    --(145.45, 571.40)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (146.90, 571.21)
    --(147.34, 569.70)
    --(153.95, 569.70)
    --(153.51, 571.28)
    --(146.90, 571.34)
    --(146.90, 571.21)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (154.71, 571.40)
    --(155.46, 568.38)
    --(160.31, 568.57)
    --(159.62, 571.72)
    --(154.71, 571.72)
    --(154.71, 571.40)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (153.70, 575.37)
    --(154.52, 571.91)
    --(159.50, 572.03)
    --(158.49, 576.00)
    --(153.64, 575.87)
    --(153.70, 575.37)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (152.76, 579.02)
    --(153.32, 576.13)
    --(158.43, 576.19)
    --(157.61, 579.28)
    --(152.69, 579.40)
    --(152.76, 579.02)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (148.60, 575.75)
    --(153.13, 575.75)
    --(153.89, 572.28)
    --(149.35, 572.28)
    --(148.60, 575.75)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (143.50, 575.87)
    --(147.91, 575.75)
    --(148.85, 572.54)
    --(144.19, 572.41)
    --(143.43, 575.69)
    --(143.50, 575.87)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (155.78, 570.96)
    --(156.16, 569.45)
    --(159.62, 569.45)
    --(159.24, 571.09)
    --(155.72, 571.02)
    --(155.78, 570.96)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (155.15, 574.93)
    --(155.40, 573.35)
    --(158.74, 573.35)
    --(158.30, 575.31)
    --(155.02, 575.18)
    --(155.15, 574.93)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (154.27, 578.65)
    --(154.46, 577.26)
    --(157.73, 577.26)
    --(157.42, 578.77)
    --(154.27, 578.77)
    --(154.27, 578.65)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (133.48, 563.84)
    --(134.43, 565.10);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (132.85, 567.69)
    --(134.17, 567.06);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (137.64, 567.75)
    --(136.94, 567.06);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (138.46, 563.91)
    --(137.26, 565.17);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (138.90, 563.65)
    --(139.84, 565.61);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (138.14, 567.69)
    --(139.46, 567.18);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (143.31, 567.69)
    --(142.49, 567.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (144.06, 563.65)
    --(142.93, 565.54);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (144.69, 563.59)
    --(145.64, 564.66);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (143.43, 567.69)
    --(145.07, 567.06);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (149.48, 567.75)
    --(149.17, 567.06);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (150.74, 563.59)
    --(149.61, 564.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (151.50, 563.78)
    --(152.44, 565.04);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (150.43, 567.81)
    --(151.94, 567.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (155.15, 567.87)
    --(154.71, 567.37);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (155.91, 563.78)
    --(155.21, 565.10);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (156.66, 563.97)
    --(157.48, 565.42);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (155.72, 567.87)
    --(157.10, 567.12);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (160.63, 567.87)
    --(160.38, 567.18);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (161.64, 564.09)
    --(160.69, 565.42);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (155.46, 568.38)
    --(156.16, 569.51);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (154.77, 571.65)
    --(155.84, 571.02);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (159.69, 571.65)
    --(159.31, 571.02);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (160.38, 568.50)
    --(159.69, 569.51);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (146.83, 559.37)
    --(148.28, 560.38);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (145.89, 563.15)
    --(147.84, 562.58);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (151.94, 559.43)
    --(151.12, 560.38);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (150.99, 563.21)
    --(150.68, 562.77);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (151.56, 563.40)
    --(153.20, 562.96);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (152.57, 559.37)
    --(153.70, 561.01);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (157.48, 559.43)
    --(156.54, 561.01);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (156.35, 563.46)
    --(155.91, 563.02);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (156.72, 563.72)
    --(158.30, 562.96);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (157.92, 559.43)
    --(158.74, 561.01);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (162.83, 559.50)
    --(161.64, 561.01);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (161.76, 563.78)
    --(161.32, 563.02);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (144.69, 575.12)
    --(147.65, 575.12)
    --(148.03, 573.23)
    --(145.01, 573.17)
    --(144.63, 575.12)
    --(144.69, 575.12)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (149.67, 575.06)
    --(152.76, 575.12)
    --(153.07, 573.04)
    --(150.24, 572.98)
    --(149.67, 574.93)
    --(149.67, 575.06)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (145.01, 579.09)
    --(145.70, 576.19)
    --(152.88, 576.19)
    --(152.57, 578.08)
    --(152.25, 579.59)
    --(146.20, 579.53)
    --(145.01, 579.21)
    --(145.01, 579.09)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (146.14, 579.46)
    --(146.46, 577.95)
    --(152.50, 578.08);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (139.84, 579.34)
    --(140.47, 575.87)
    --(145.39, 576.00)
    --(144.76, 579.09)
    --(144.25, 579.53)
    --(140.60, 579.46)
    --(139.84, 579.40)
    --(139.84, 579.34)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (140.60, 579.46)
    --(141.04, 577.51)
    --(144.38, 577.70)
    --(144.25, 579.53);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (146.58, 577.89)
    --(145.76, 576.25);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (141.04, 577.51)
    --(140.60, 576.00);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (144.50, 577.57)
    --(145.32, 576.19);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (138.08, 575.56)
    --(138.77, 572.09)
    --(143.56, 572.16)
    --(142.68, 575.69)
    --(138.14, 575.75);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (139.28, 575.31)
    --(139.53, 573.48)
    --(142.55, 573.54)
    --(142.24, 575.18)
    --(139.34, 575.18);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (139.34, 575.06)
    --(138.08, 575.69);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (139.59, 573.48)
    --(138.83, 572.03);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (142.55, 573.48)
    --(143.50, 572.22);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (142.24, 575.18)
    --(142.55, 575.69);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (140.03, 571.65)
    --(140.85, 568.00)
    --(145.89, 568.06)
    --(144.88, 571.78)
    --(140.03, 571.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (140.91, 568.00)
    --(141.73, 569.70)
    --(144.63, 569.76)
    --(145.89, 568.13);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (141.73, 569.70)
    --(141.42, 571.34)
    --(144.69, 571.46)
    --(144.69, 569.83);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (141.35, 571.34)
    --(140.22, 571.65);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (144.69, 571.34)
    --(144.88, 571.78);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (134.99, 579.34)
    --(135.56, 575.94)
    --(140.35, 575.94)
    --(139.72, 579.59)
    --(134.93, 579.46);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (135.62, 576.00)
    --(136.50, 577.20)
    --(139.53, 577.20)
    --(140.22, 576.00);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (136.50, 577.20)
    --(136.06, 579.21)
    --(134.99, 579.46);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (136.13, 579.15)
    --(139.34, 579.21)
    --(139.78, 579.46);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (139.34, 579.09)
    --(139.46, 577.20);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (133.10, 575.62)
    --(133.61, 571.84)
    --(138.77, 571.91)
    --(137.83, 575.69)
    --(133.04, 575.75);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (133.67, 571.91)
    --(134.55, 573.61)
    --(137.57, 573.67)
    --(138.71, 571.97);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (134.61, 573.67)
    --(134.24, 575.18)
    --(133.17, 575.69);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (134.36, 575.24)
    --(137.64, 575.24)
    --(137.89, 575.69);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (137.64, 575.18)
    --(137.64, 573.67);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (134.87, 571.59)
    --(135.50, 568.00)
    --(140.66, 567.94)
    --(139.72, 571.72)
    --(134.87, 571.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (135.50, 567.94)
    --(136.50, 569.51)
    --(139.34, 569.64)
    --(140.72, 567.94);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (136.57, 569.64)
    --(136.50, 570.83)
    --(134.99, 571.65);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (136.57, 570.90)
    --(139.34, 570.90)
    --(139.72, 571.65);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (139.46, 570.83)
    --(139.46, 569.64);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (127.37, 567.37)
    --(127.69, 563.53)
    --(133.17, 563.59)
    --(132.54, 567.43)
    --(127.31, 567.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (127.50, 566.93);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (127.69, 563.59)
    --(129.01, 565.17)
    --(131.91, 565.17)
    --(133.17, 563.46);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (129.13, 565.17)
    --(128.82, 566.49)
    --(127.43, 567.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (128.82, 566.61)
    --(131.72, 566.68)
    --(132.54, 567.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (131.97, 565.17)
    --(131.72, 566.68);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (131.72, 566.49);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (129.51, 571.59)
    --(129.83, 567.75)
    --(134.99, 568.00)
    --(134.49, 571.84)
    --(129.51, 571.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (129.89, 567.87)
    --(131.21, 569.13)
    --(133.98, 569.13)
    --(134.99, 568.00)
    --(134.87, 568.00);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (129.51, 571.72)
    --(131.15, 570.90)
    --(131.21, 569.13);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (131.21, 570.90)
    --(133.86, 571.02)
    --(134.49, 571.78);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (133.98, 570.96)
    --(133.86, 569.20);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (129.83, 579.59)
    --(130.33, 575.94)
    --(135.12, 576.00)
    --(134.68, 579.59)
    --(129.76, 579.59);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (130.46, 575.87)
    --(131.34, 577.57)
    --(134.11, 577.64)
    --(135.06, 575.94);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (129.83, 579.46)
    --(131.28, 578.71)
    --(131.40, 577.57);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (131.28, 578.71)
    --(134.17, 578.71)
    --(134.17, 577.57);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (134.17, 578.71)
    --(134.74, 579.59);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (127.69, 575.56)
    --(128.19, 571.78)
    --(133.48, 571.91)
    --(132.91, 575.69)
    --(127.69, 575.69);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (128.25, 571.84)
    --(129.32, 573.48)
    --(132.35, 573.48)
    --(133.54, 571.91);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (129.39, 573.48)
    --(129.32, 574.74)
    --(127.69, 575.69);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (129.32, 574.80)
    --(132.22, 574.87)
    --(132.98, 575.69);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (132.22, 574.80)
    --(132.22, 573.48);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (146.46, 568.25)
    --(147.46, 569.70);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (144.19, 572.35)
    --(145.01, 573.17);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (145.39, 571.84)
    --(146.96, 571.28);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (154.02, 571.84)
    --(153.45, 571.21);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (148.85, 572.47)
    --(148.03, 573.23);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (147.97, 575.81)
    --(147.72, 575.06);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (143.50, 575.75)
    --(144.69, 575.06);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (148.66, 575.81)
    --(149.67, 575.06);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (149.35, 572.28)
    --(150.24, 572.98);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (154.02, 572.28)
    --(153.07, 572.98);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (153.13, 575.81)
    --(152.82, 575.06);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (152.94, 575.06);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (153.39, 576.13)
    --(154.46, 577.26);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (152.82, 579.34)
    --(154.39, 578.71);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (157.67, 579.28)
    --(157.42, 578.77);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (158.49, 576.13)
    --(157.73, 577.26);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (154.90, 568.44)
    --(153.95, 569.70);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (154.71, 571.91)
    --(155.53, 573.35);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (153.76, 575.81)
    --(155.21, 575.12);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (158.61, 575.94)
    --(158.36, 575.18);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (159.56, 572.03)
    --(158.74, 573.35);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (124.28, 571.59)
    --(124.54, 567.87)
    --(129.76, 567.75)
    --(129.26, 571.72)
    --(124.28, 571.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (124.66, 567.87)
    --(125.54, 569.13)
    --(128.69, 569.20)
    --(129.70, 567.81);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (125.54, 569.13)
    --(125.54, 570.90)
    --(124.47, 571.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (125.54, 570.83)
    --(128.63, 570.83)
    --(129.20, 571.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (128.63, 570.83)
    --(128.69, 569.13);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (121.76, 567.56)
    --(122.14, 563.34)
    --(127.50, 563.46)
    --(127.18, 567.50)
    --(121.76, 567.56);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (122.08, 563.40)
    --(123.15, 564.72)
    --(126.36, 564.66)
    --(127.50, 563.46);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (126.43, 564.72)
    --(126.36, 566.61)
    --(127.18, 567.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (126.36, 566.68)
    --(122.96, 566.61)
    --(123.09, 564.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (123.02, 566.55)
    --(121.83, 567.56);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (118.93, 571.72)
    --(119.12, 567.81)
    --(124.41, 567.81)
    --(123.97, 571.72)
    --(118.93, 571.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (119.12, 567.87)
    --(120.19, 569.13)
    --(123.09, 569.20)
    --(124.47, 567.87);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (120.25, 569.13)
    --(120.19, 570.90)
    --(118.93, 571.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (120.31, 570.90)
    --(123.09, 570.96)
    --(123.97, 571.65);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (123.09, 570.90)
    --(123.09, 569.13);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (116.41, 567.56)
    --(116.41, 563.53)
    --(121.89, 563.53)
    --(121.57, 567.50)
    --(116.41, 567.56);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (116.41, 563.53)
    --(117.61, 564.91)
    --(120.82, 564.91)
    --(121.76, 563.53);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (117.61, 564.98)
    --(117.48, 566.49)
    --(116.41, 567.56);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (117.42, 566.61)
    --(120.88, 566.55)
    --(121.57, 567.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (120.88, 566.55)
    --(120.82, 564.91);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (122.46, 575.62)
    --(122.90, 571.91)
    --(127.94, 571.91)
    --(127.37, 575.69)
    --(122.46, 575.62);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (123.02, 571.91)
    --(123.84, 573.29)
    --(126.93, 573.29)
    --(127.94, 571.97);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (123.84, 573.35)
    --(123.72, 574.87)
    --(122.58, 575.62);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (123.78, 574.93)
    --(127.06, 574.99)
    --(127.37, 575.69);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (127.06, 574.87)
    --(126.93, 573.35);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (124.66, 579.40)
    --(125.04, 576.00)
    --(129.89, 576.00)
    --(129.39, 579.59)
    --(124.66, 579.59)
    --(124.66, 579.40);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (125.17, 576.00)
    --(125.92, 577.20)
    --(129.01, 577.20)
    --(129.95, 576.06);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (129.07, 577.20)
    --(128.88, 578.83)
    --(129.32, 579.59);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (128.82, 578.83)
    --(125.61, 578.83)
    --(125.92, 577.20);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (125.54, 578.90)
    --(124.66, 579.59);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (119.50, 579.84)
    --(119.81, 576.06)
    --(124.66, 576.06)
    --(124.09, 579.78)
    --(119.50, 579.78);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (119.81, 576.06)
    --(120.82, 577.20)
    --(123.59, 577.26)
    --(124.66, 576.13);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (119.50, 579.72)
    --(120.57, 579.21)
    --(120.82, 577.13);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (120.69, 579.21)
    --(123.59, 579.21)
    --(124.16, 579.84);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (123.72, 579.21)
    --(123.65, 577.20)
    --(123.65, 577.26);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (117.10, 575.69)
    --(117.35, 571.91)
    --(122.71, 571.91)
    --(122.08, 575.69)
    --(117.17, 575.62);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (117.35, 571.91)
    --(118.43, 573.23)
    --(121.57, 573.23)
    --(122.71, 572.03);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (118.49, 573.23)
    --(118.36, 574.68)
    --(117.23, 575.62);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (118.43, 574.68)
    --(121.57, 574.68)
    --(122.08, 575.69);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (121.57, 574.61)
    --(121.51, 573.23);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (114.39, 579.46)
    --(114.46, 576.06)
    --(119.50, 576.19)
    --(119.12, 579.65)
    --(114.39, 579.59);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (114.46, 576.13)
    --(115.28, 577.13)
    --(118.68, 577.20)
    --(119.43, 576.19);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (115.34, 577.07)
    --(115.21, 579.02)
    --(114.39, 579.65);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (115.34, 579.02)
    --(118.68, 579.09)
    --(119.12, 579.59);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (118.68, 578.96)
    --(118.74, 577.20);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (111.94, 575.62)
    --(112.06, 571.97)
    --(117.04, 571.97)
    --(116.72, 575.69)
    --(111.94, 575.75);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (112.06, 571.97)
    --(112.94, 572.91)
    --(116.22, 572.98)
    --(117.04, 571.97);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (112.94, 572.98)
    --(112.76, 574.87)
    --(112.06, 575.69);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (112.76, 574.87)
    --(116.22, 574.87)
    --(116.72, 575.69);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (116.22, 574.87)
    --(116.22, 572.98);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (113.26, 571.65)
    --(113.20, 568.13)
    --(118.80, 568.00)
    --(118.49, 571.72)
    --(113.26, 571.78);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (113.26, 568.00)
    --(114.27, 568.94)
    --(117.92, 568.94)
    --(118.74, 568.06);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (114.27, 568.94)
    --(114.27, 570.90)
    --(113.26, 571.65);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (114.08, 571.02)
    --(117.98, 571.02)
    --(118.43, 571.65);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (118.49, 571.59);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (118.05, 570.90)
    --(117.98, 568.94);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (110.55, 567.56)
    --(110.80, 563.53)
    --(116.03, 563.59)
    --(115.72, 567.69)
    --(110.55, 567.69);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (110.87, 563.59)
    --(111.81, 564.66)
    --(115.21, 564.66)
    --(115.91, 563.65);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (111.87, 564.66)
    --(111.75, 566.74)
    --(110.55, 567.69);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (111.81, 566.74)
    --(115.09, 566.74)
    --(115.72, 567.62);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (115.09, 566.74)
    --(115.15, 564.66);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (109.10, 579.46)
    --(109.35, 576.25)
    --(114.14, 576.25)
    --(113.70, 579.59)
    --(109.04, 579.59);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (109.35, 576.38)
    --(110.17, 577.20)
    --(113.32, 577.20)
    --(114.08, 576.31);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (110.17, 577.20)
    --(110.11, 578.77)
    --(109.10, 579.53);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (110.17, 578.83)
    --(113.26, 578.83)
    --(113.64, 579.53);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (113.32, 578.77)
    --(113.32, 577.20);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (106.52, 575.75)
    --(106.65, 572.09)
    --(111.81, 571.97)
    --(111.43, 575.69)
    --(106.46, 575.75);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (106.71, 572.09)
    --(107.59, 573.35)
    --(110.74, 573.42)
    --(111.81, 572.03);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (107.72, 573.48)
    --(107.65, 574.87)
    --(106.52, 575.75);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (107.72, 574.87)
    --(110.80, 574.87)
    --(111.43, 575.69);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (110.80, 574.80)
    --(110.80, 573.42);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (107.78, 571.65)
    --(107.91, 567.94)
    --(112.94, 568.00)
    --(112.82, 571.59)
    --(107.78, 571.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (107.91, 567.87)
    --(108.85, 568.88)
    --(112.25, 568.88)
    --(112.94, 568.00);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (108.91, 568.94)
    --(108.91, 570.90)
    --(107.78, 571.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (108.98, 570.90)
    --(112.38, 570.90)
    --(112.88, 571.59);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (112.38, 570.90)
    --(112.25, 568.88);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (105.13, 567.56)
    --(105.26, 563.46)
    --(110.49, 563.53)
    --(110.30, 567.56)
    --(105.13, 567.56);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (105.26, 563.46)
    --(106.39, 564.72)
    --(109.61, 564.79)
    --(110.49, 563.53);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (106.52, 564.72)
    --(106.39, 566.68)
    --(105.13, 567.56);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (106.39, 566.74)
    --(109.61, 566.74)
    --(110.30, 567.62);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (109.61, 566.68)
    --(109.61, 564.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (104.00, 579.72)
    --(104.00, 576.00)
    --(108.72, 576.00)
    --(108.41, 579.78)
    --(104.06, 579.78);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (103.94, 576.06)
    --(104.88, 577.13)
    --(107.78, 577.13)
    --(108.72, 576.06);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (105.01, 577.13)
    --(105.01, 579.09)
    --(104.00, 579.84);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (105.01, 579.15)
    --(107.84, 579.15)
    --(108.47, 579.78);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (107.84, 579.15)
    --(107.78, 577.07);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (101.35, 575.75)
    --(101.29, 572.03)
    --(106.27, 572.03)
    --(106.14, 575.75)
    --(101.35, 575.81);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (101.29, 572.03)
    --(102.30, 573.23)
    --(105.39, 573.23)
    --(106.27, 572.09);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (102.36, 573.23)
    --(102.43, 574.99)
    --(101.35, 575.75);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (102.55, 575.12)
    --(105.45, 575.12)
    --(106.14, 575.75);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (105.51, 574.99)
    --(105.45, 573.23);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (102.43, 571.78)
    --(102.43, 567.87)
    --(107.40, 567.87)
    --(107.40, 571.72)
    --(102.43, 571.78);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (102.43, 567.94)
    --(103.31, 568.94)
    --(106.71, 568.94)
    --(107.40, 567.94);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (103.43, 568.94)
    --(103.37, 570.77)
    --(102.43, 571.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (103.43, 570.77)
    --(106.71, 570.77)
    --(107.40, 571.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (106.71, 570.77)
    --(106.71, 568.94);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (99.65, 567.56)
    --(99.65, 563.46)
    --(104.82, 563.46)
    --(104.82, 567.56)
    --(99.65, 567.56)
    --(99.72, 567.62);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (99.72, 563.59)
    --(100.72, 565.17)
    --(104.00, 565.17)
    --(104.82, 563.46);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (100.79, 565.17)
    --(100.85, 566.74)
    --(99.72, 567.56);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (100.85, 566.74)
    --(103.81, 566.74)
    --(104.76, 567.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (103.87, 566.74)
    --(103.94, 565.17);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (98.27, 576.06)
    --(103.62, 576.06)
    --(103.43, 579.65)
    --(98.39, 579.72)
    --(98.46, 576.19);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (98.46, 576.06)
    --(99.72, 577.45)
    --(102.55, 577.51)
    --(103.62, 576.19);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (99.72, 577.45)
    --(99.65, 579.46)
    --(98.39, 579.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (99.72, 579.40)
    --(102.87, 579.34)
    --(103.50, 579.65);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (103.50, 579.46);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (102.93, 579.28)
    --(102.55, 577.45);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (93.92, 567.69)
    --(94.05, 563.46)
    --(99.21, 563.46)
    --(99.21, 567.69)
    --(93.92, 567.81);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (94.11, 563.46)
    --(95.12, 564.72)
    --(98.08, 564.72)
    --(99.28, 563.53);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (98.08, 564.72)
    --(98.08, 567.12)
    --(99.28, 567.69);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (97.95, 567.18)
    --(94.99, 567.12)
    --(94.05, 567.69);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (95.06, 567.06)
    --(95.24, 564.66);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (96.88, 571.59)
    --(96.94, 568.00)
    --(102.17, 568.00)
    --(101.86, 571.72)
    --(96.88, 571.78);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (97.01, 568.00)
    --(97.83, 569.13)
    --(101.29, 569.13)
    --(102.17, 568.00);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (101.29, 569.20)
    --(101.29, 571.46)
    --(101.86, 571.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (97.83, 569.20)
    --(97.76, 571.46)
    --(96.94, 571.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (97.89, 571.40)
    --(101.29, 571.40);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (95.94, 575.50)
    --(95.94, 572.22)
    --(100.98, 572.22)
    --(100.98, 575.69)
    --(95.87, 575.69);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (95.87, 572.22)
    --(96.63, 573.42)
    --(100.16, 573.42)
    --(101.04, 572.22);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (100.16, 573.48)
    --(100.16, 575.31)
    --(100.98, 575.69);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (100.16, 575.31)
    --(96.69, 575.37)
    --(95.81, 575.69);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (96.63, 573.48)
    --(96.63, 575.37);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (93.29, 579.84)
    --(93.29, 576.06)
    --(98.27, 576.00)
    --(98.08, 579.78)
    --(93.42, 579.78);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (93.29, 576.06)
    --(94.30, 577.45)
    --(97.20, 577.45)
    --(98.20, 576.06);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (97.20, 577.45)
    --(97.26, 579.53)
    --(98.14, 579.78);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (97.26, 579.53)
    --(94.36, 579.53)
    --(93.29, 579.78);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (94.24, 579.53)
    --(94.36, 577.45);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (88.31, 579.72)
    --(88.19, 576.13)
    --(92.91, 576.13)
    --(92.85, 579.78)
    --(88.31, 579.78);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (88.31, 576.19)
    --(89.07, 577.45)
    --(92.09, 577.51)
    --(92.91, 576.13);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (92.09, 577.51)
    --(92.09, 579.53)
    --(92.91, 579.78);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (92.09, 579.53)
    --(89.20, 579.53)
    --(89.20, 577.51);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (89.20, 579.53)
    --(88.25, 579.78);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (90.52, 575.69)
    --(90.39, 571.84)
    --(95.43, 571.84)
    --(95.56, 575.75)
    --(90.52, 575.87);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (90.33, 571.91)
    --(91.40, 573.29)
    --(94.55, 573.29)
    --(95.37, 571.84);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (94.61, 573.29)
    --(94.61, 575.24)
    --(95.56, 575.75);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (94.55, 575.24)
    --(91.28, 575.24)
    --(91.34, 573.35);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (91.34, 575.18)
    --(90.52, 575.81);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (91.59, 571.59)
    --(91.59, 568.00)
    --(96.63, 568.00)
    --(96.63, 571.72)
    --(91.59, 571.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (91.65, 568.06)
    --(92.66, 569.26)
    --(95.87, 569.26)
    --(96.57, 568.13);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (95.94, 569.26)
    --(95.87, 571.40)
    --(96.57, 571.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (96.00, 571.40)
    --(92.41, 571.46)
    --(92.66, 569.20);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (92.47, 571.46)
    --(91.59, 571.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (88.57, 567.31)
    --(88.57, 563.65)
    --(93.67, 563.65)
    --(93.67, 567.50)
    --(88.57, 567.62)
    --(88.57, 567.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (88.57, 566.68);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (88.57, 563.65)
    --(89.51, 564.98)
    --(92.79, 564.98)
    --(93.61, 563.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (93.42, 563.97);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (92.79, 564.91)
    --(92.91, 567.12)
    --(93.73, 567.56);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (92.91, 567.12)
    --(89.32, 567.12)
    --(89.45, 565.04);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (89.45, 567.06)
    --(88.57, 567.62);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (85.92, 571.53)
    --(85.98, 567.87)
    --(91.40, 567.81)
    --(91.15, 571.65)
    --(85.92, 571.65);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (86.05, 567.81)
    --(86.93, 569.26)
    --(90.39, 569.26)
    --(91.40, 567.87);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (90.46, 569.26)
    --(90.39, 571.21)
    --(91.21, 571.65);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (90.46, 571.21)
    --(86.74, 571.15)
    --(86.93, 569.32);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (86.80, 571.15)
    --(85.92, 571.59);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (85.10, 575.50)
    --(85.04, 571.78)
    --(90.33, 571.84)
    --(90.27, 575.81)
    --(85.10, 575.81)
    --(85.10, 575.50);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (85.17, 571.84)
    --(85.92, 573.23)
    --(89.13, 573.23)
    --(90.39, 571.84);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (89.20, 573.23)
    --(89.13, 575.24)
    --(90.33, 575.75);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (89.20, 575.18)
    --(85.92, 575.18)
    --(85.98, 573.23);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (86.05, 575.24)
    --(85.10, 575.87);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (82.96, 579.72)
    --(82.77, 576.19)
    --(87.81, 576.13)
    --(87.69, 579.72)
    --(82.96, 579.78);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (82.83, 576.25)
    --(83.84, 577.45)
    --(86.61, 577.45)
    --(87.81, 576.19);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (86.55, 577.45)
    --(86.61, 579.21)
    --(87.75, 579.78);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (86.61, 579.21)
    --(83.97, 579.21)
    --(83.91, 577.45);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (83.97, 579.15)
    --(83.02, 579.72);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (83.02, 579.84);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (78.87, 567.31)
    --(78.30, 563.46)
    --(88.06, 563.65)
    --(88.25, 567.43)
    --(78.80, 567.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (78.36, 563.53)
    --(79.37, 564.91)
    --(86.43, 565.04)
    --(88.13, 563.65);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (86.49, 564.98)
    --(86.74, 566.80)
    --(88.19, 567.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (88.25, 567.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (86.74, 566.80)
    --(79.50, 566.87)
    --(79.37, 564.91);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (79.50, 566.80)
    --(78.93, 567.37);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (79.18, 571.59)
    --(78.93, 567.75)
    --(85.73, 567.75)
    --(85.54, 571.59)
    --(79.18, 571.59);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (79.12, 567.81)
    --(80.19, 568.88)
    --(84.72, 568.82)
    --(85.73, 567.87);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (84.79, 568.76)
    --(84.79, 570.90)
    --(85.54, 571.59);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (84.72, 570.83)
    --(80.13, 570.83)
    --(80.19, 568.76);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (80.19, 570.83)
    --(79.31, 571.59);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (79.43, 571.97)
    --(84.91, 571.97)
    --(84.72, 575.75)
    --(80.13, 575.81)
    --(79.31, 571.91);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (79.50, 571.97)
    --(80.50, 573.67)
    --(82.90, 573.73)
    --(84.98, 571.97);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (80.44, 573.67)
    --(80.82, 575.31)
    --(80.13, 575.81);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (80.82, 575.31)
    --(82.90, 575.31)
    --(82.90, 573.73);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (82.90, 575.37)
    --(84.60, 575.69);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (80.31, 580.09)
    --(82.83, 580.09)
    --(82.58, 575.94)
    --(80.06, 575.94)
    --(80.38, 580.09);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (80.82, 582.80)
    --(80.50, 580.35)
    --(85.48, 580.22)
    --(85.86, 582.68)
    --(80.82, 582.68);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (80.57, 580.47)
    --(81.13, 581.61)
    --(84.72, 581.61)
    --(85.48, 580.28);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (84.79, 581.61)
    --(84.85, 582.61);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (81.20, 581.67)
    --(81.51, 582.68);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (86.68, 582.61)
    --(86.49, 580.22)
    --(91.65, 580.22)
    --(91.97, 582.49)
    --(86.80, 582.55);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (86.43, 580.16)
    --(87.24, 580.91)
    --(91.09, 581.04)
    --(91.78, 580.22);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (87.24, 581.04)
    --(87.24, 582.43)
    --(87.31, 582.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (91.09, 580.91)
    --(91.21, 582.49);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (92.47, 582.36)
    --(92.41, 580.28)
    --(96.44, 580.22)
    --(96.63, 582.49)
    --(92.41, 582.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (92.41, 580.22)
    --(93.23, 581.10)
    --(95.56, 581.10)
    --(96.44, 580.22);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (95.62, 581.17)
    --(95.69, 582.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (93.23, 581.17)
    --(93.29, 582.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (97.07, 582.36)
    --(97.07, 580.03)
    --(101.29, 580.03)
    --(101.35, 582.36)
    --(97.13, 582.49);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (97.20, 580.09)
    --(97.70, 580.91)
    --(100.54, 580.91)
    --(101.42, 580.03);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (97.76, 580.91)
    --(97.83, 582.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (100.54, 580.91)
    --(100.54, 582.30);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (101.86, 582.36)
    --(101.80, 580.03)
    --(106.08, 580.03)
    --(106.08, 582.49)
    --(101.86, 582.55);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (101.98, 582.55);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (101.86, 580.03)
    --(102.74, 581.04)
    --(105.39, 580.98)
    --(106.08, 580.09);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (102.74, 580.98)
    --(102.68, 582.55);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (105.32, 580.98)
    --(105.39, 582.49);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (107.02, 582.43)
    --(106.96, 580.03)
    --(111.06, 580.03)
    --(111.06, 582.49)
    --(106.96, 582.49);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (106.96, 580.03)
    --(107.53, 581.04)
    --(110.24, 581.04)
    --(111.12, 580.09);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (107.53, 581.04)
    --(107.53, 582.49);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (110.17, 581.04)
    --(110.17, 582.49);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (111.37, 582.43)
    --(111.37, 579.97)
    --(115.59, 579.97)
    --(115.59, 582.43)
    --(111.43, 582.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (111.37, 582.55);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (111.37, 579.91)
    --(111.87, 580.72)
    --(114.90, 580.72)
    --(115.59, 579.97);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (114.90, 580.79)
    --(114.90, 582.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (111.87, 580.72)
    --(111.87, 582.49);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (115.97, 582.36)
    --(115.97, 579.91)
    --(119.94, 579.97)
    --(119.87, 582.49)
    --(115.91, 582.49);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (115.97, 579.91)
    --(116.54, 580.72)
    --(119.18, 580.72)
    --(119.94, 579.97);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (116.60, 580.72)
    --(116.60, 582.49);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (119.18, 580.72)
    --(119.18, 582.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (120.38, 582.43)
    --(120.31, 580.03)
    --(124.60, 579.97)
    --(124.60, 582.43)
    --(120.38, 582.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (120.44, 579.97)
    --(121.07, 580.72)
    --(124.09, 580.72)
    --(124.54, 579.97);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (125.10, 582.43)
    --(125.10, 579.91)
    --(129.51, 579.91)
    --(129.51, 582.49)
    --(125.10, 582.49);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (125.10, 579.91)
    --(125.61, 580.72)
    --(129.13, 580.72)
    --(129.51, 579.97);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (125.67, 580.72)
    --(125.61, 582.55);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (129.13, 580.72)
    --(129.13, 582.55);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (121.07, 580.72)
    --(121.07, 582.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (124.22, 580.72)
    --(124.22, 582.43);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (129.95, 582.49)
    --(129.95, 579.91)
    --(134.43, 579.97)
    --(134.43, 582.49)
    --(129.95, 582.61);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (129.95, 579.97)
    --(130.65, 580.72)
    --(133.73, 580.72)
    --(134.43, 580.03);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (130.77, 580.72)
    --(130.77, 582.49);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (133.80, 580.72)
    --(133.80, 582.55);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (134.74, 582.43)
    --(134.74, 579.91)
    --(139.40, 579.91)
    --(139.28, 582.55)
    --(134.68, 582.49);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (134.80, 579.91)
    --(135.37, 580.60)
    --(138.65, 580.60)
    --(139.46, 579.91);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (135.50, 580.66)
    --(135.50, 582.49);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (138.65, 580.60)
    --(138.65, 582.55);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (139.65, 582.61)
    --(139.65, 579.91)
    --(143.87, 579.91)
    --(143.75, 582.61)
    --(139.65, 582.68);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (139.72, 579.91)
    --(140.47, 580.72)
    --(143.18, 580.72)
    --(143.87, 579.91);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (143.18, 580.72)
    --(143.18, 582.61);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (140.54, 580.72)
    --(140.47, 582.61);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (144.31, 582.61)
    --(144.31, 579.91)
    --(148.54, 579.91)
    --(148.54, 582.61)
    --(144.31, 582.68);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (144.31, 579.84)
    --(145.13, 580.79)
    --(147.91, 580.79)
    --(148.54, 579.91);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (145.20, 580.72)
    --(145.20, 582.61);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (147.91, 580.72)
    --(147.91, 582.61);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (148.85, 582.55)
    --(148.79, 579.84)
    --(152.57, 579.84)
    --(152.57, 582.49)
    --(148.85, 582.49);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (148.85, 579.91)
    --(149.54, 580.79)
    --(152.13, 580.79)
    --(152.50, 579.91);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (149.67, 580.79)
    --(149.67, 582.49);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (152.25, 580.72)
    --(152.25, 582.49);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (153.01, 582.36)
    --(153.01, 579.78)
    --(157.61, 579.78)
    --(156.98, 582.55)
    --(153.01, 582.49);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (153.13, 579.78)
    --(153.83, 580.66)
    --(157.04, 580.66)
    --(157.61, 579.84);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (157.10, 580.66)
    --(156.60, 582.49);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (153.95, 580.72)
    --(153.76, 582.49);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (155.53, 584.82)
    --(157.17, 584.82)
    --(157.17, 584.00)
    --(155.53, 584.00)
    --(155.53, 584.82)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (107.21, 648.25)
    --(107.46, 646.30)
    --(130.52, 646.49)
    --(130.58, 648.38);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (125.04, 562.77)
    --(126.49, 563.34);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin ellipse object
  pickup pencircle scaled 0.94;
  path c;
  c = fullcircle scaled 22.80 yscaled 1.00
         rotated 0.00 shifted (458.33,634.27);
  draw c withcolor (0.00,0.00,0.00);
% End ellipse object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt 2 etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((454.55,628.60))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt Fin3D/Troot etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((721.89,494.36))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt Truss/Deflection etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((721.89,437.67))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt Fin2D/Troot etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((721.89,466.02))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt Truss/Stress etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((721.89,409.32))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (727.56, 522.71)
    --(713.39, 522.71)
    --(713.39, 395.15)
    --(727.56, 395.15);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt Fin3D/Troot etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((719.06,338.46))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt Truss/Deflection etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((719.06,281.76))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt Fin2D/Troot etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((719.06,310.11))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt Truss/Stress etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((719.06,253.42))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (727.56, 366.80)
    --(713.39, 366.80)
    --(713.39, 239.24)
    --(727.56, 239.24);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt Fin3D/Troot etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((721.89,635.15))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt Truss/Deflection etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((721.89,578.46))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt Fin2D/Troot etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((721.89,606.80))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt Truss/Stress etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((721.89,550.11))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (727.56, 663.50)
    --(713.39, 663.50)
    --(713.39, 535.94)
    --(727.56, 535.94);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin ellipse object
  pickup pencircle scaled 0.94;
  path c;
  c = fullcircle scaled 22.80 yscaled 1.00
         rotated 0.00 shifted (466.71,319.50);
  draw c withcolor (0.00,0.00,0.00);
% End ellipse object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt 2 etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((462.93,313.83))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin ellipse object
  pickup pencircle scaled 0.94;
  path c;
  c = fullcircle scaled 22.80 yscaled 1.00
         rotated 0.00 shifted (498.02,483.91);
  draw c withcolor (0.00,0.00,0.00);
% End ellipse object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt 2 etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((494.24,478.24))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (258.46, 699.91)
    --(260.72, 699.91);
  draw p withcolor (0.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (265.83, 699.91)
    --(268.09, 699.91);
  draw p withcolor (1.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (269.54, 699.91)
    --(271.81, 699.91);
  draw p withcolor (0.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (276.47, 699.91)
    --(278.74, 699.91);
  draw p withcolor (0.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (283.84, 699.91)
    --(286.05, 699.91);
  draw p withcolor (1.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (287.56, 699.91)
    --(289.83, 699.91);
  draw p withcolor (0.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (294.49, 699.91)
    --(296.76, 699.91);
  draw p withcolor (0.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (301.86, 699.91)
    --(304.06, 699.91);
  draw p withcolor (1.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (305.57, 699.91)
    --(307.78, 699.91);
  draw p withcolor (0.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (258.46, 693.17)
    --(260.72, 693.17);
  draw p withcolor (0.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (265.83, 693.17)
    --(268.09, 693.17);
  draw p withcolor (1.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (269.54, 693.17)
    --(271.81, 693.17);
  draw p withcolor (0.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (276.47, 693.17)
    --(278.74, 693.17);
  draw p withcolor (0.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (283.84, 693.17)
    --(286.05, 693.17);
  draw p withcolor (1.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (287.56, 693.17)
    --(289.83, 693.17);
  draw p withcolor (0.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (294.49, 693.17)
    --(296.76, 693.17);
  draw p withcolor (0.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (301.86, 693.17)
    --(304.06, 693.17);
  draw p withcolor (1.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (305.57, 693.17)
    --(307.78, 693.17);
  draw p withcolor (0.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (246.80, 696.13)
    --(249.07, 696.13);
  draw p withcolor (1.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (251.02, 696.13)
    --(253.23, 696.13);
  draw p withcolor (0.00,1.00,0.00);
% End polyline object
% Begin ellipse object
  pickup pencircle scaled 0.94;
  path c;
  c = fullcircle scaled 22.80 yscaled 1.00
         rotated 0.00 shifted (173.92,515.09);
  draw c withcolor (0.00,0.00,0.00);
% End ellipse object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt 1 etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((170.14,509.42))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin ellipse object
  pickup pencircle scaled 0.94;
  path c;
  c = fullcircle scaled 22.80 yscaled 1.00
         rotated 0.00 shifted (188.09,653.04);
  draw c withcolor (0.00,0.00,0.00);
% End ellipse object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt 1 etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((184.31,647.37))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin ellipse object
  pickup pencircle scaled 0.94;
  path c;
  c = fullcircle scaled 22.80 yscaled 1.00
         rotated 0.00 shifted (180.41,368.76);
  draw c withcolor (0.00,0.00,0.00);
% End ellipse object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt 1 etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((176.63,363.09))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin ellipse object
  pickup pencircle scaled 0.94;
  path c;
  c = fullcircle scaled 33.64 yscaled 1.00
         rotated 0.00 shifted (605.67,303.50);
  draw c withcolor (0.00,0.00,0.00);
% End ellipse object
% Begin ellipse object
  pickup pencircle scaled 0.94;
  path c;
  c = fullcircle scaled 33.64 yscaled 1.00
         rotated 0.00 shifted (604.72,456.57);
  draw c withcolor (0.00,0.00,0.00);
% End ellipse object
% Begin ellipse object
  pickup pencircle scaled 0.94;
  path c;
  c = fullcircle scaled 33.64 yscaled 1.00
         rotated 0.00 shifted (606.61,597.35);
  draw c withcolor (0.00,0.00,0.00);
% End ellipse object
% Begin ellipse object
  pickup pencircle scaled 0.94;
  path c;
  c = fullcircle scaled 33.64 yscaled 1.00
         rotated 0.00 shifted (38.74,624.76);
  draw c withcolor (0.00,0.00,0.00);
% End ellipse object
% Begin ellipse object
  pickup pencircle scaled 0.94;
  path c;
  c = fullcircle scaled 33.64 yscaled 1.00
         rotated 0.00 shifted (39.69,468.85);
  draw c withcolor (0.00,0.00,0.00);
% End ellipse object
% Begin ellipse object
  pickup pencircle scaled 0.94;
  path c;
  c = fullcircle scaled 33.64 yscaled 1.00
         rotated 0.00 shifted (39.69,324.28);
  draw c withcolor (0.00,0.00,0.00);
% End ellipse object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 1.89;
  path p;
  p = (288.19, 282.71)
    --(443.15, 282.71);
  draw p withcolor (0.00,0.00,0.00);
% Draw arrowhead type 0
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path arr;
  arr = (297.83, 284.60)
      --(290.27, 282.71)
      --(297.83, 280.82);
  draw arr withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 1.89;
  path p;
  p = (291.97, 312.00)
    --(442.20, 312.00);
  draw p withcolor (0.00,0.00,0.00);
% Draw arrowhead type 0
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path arr;
  arr = (432.57, 310.11)
      --(440.13, 312.00)
      --(432.57, 313.89);
  draw arr withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 1.89;
  path p;
  p = (278.74, 466.96)
    --(428.98, 466.96);
  draw p withcolor (0.00,0.00,0.00);
% Draw arrowhead type 0
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path arr;
  arr = (419.34, 465.07)
      --(426.90, 466.96)
      --(419.34, 468.85);
  draw arr withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 1.89;
  path p;
  p = (274.02, 593.57)
    --(428.98, 593.57);
  draw p withcolor (0.00,0.00,0.00);
% Draw arrowhead type 0
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path arr;
  arr = (283.65, 595.46)
      --(276.09, 593.57)
      --(283.65, 591.69);
  draw arr withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 1.89;
  path p;
  p = (277.80, 622.87)
    --(428.03, 622.87);
  draw p withcolor (0.00,0.00,0.00);
% Draw arrowhead type 0
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path arr;
  arr = (418.39, 620.98)
      --(425.95, 622.87)
      --(418.39, 624.76);
  draw arr withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 1.89;
  path p;
  p = (184.25, 296.88)
    --(524.41, 296.88);
  draw p withcolor (0.00,0.00,0.00) dashed evenly scaled 1.82;
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 1.89;
  path p;
  p = (171.02, 451.84)
    --(511.18, 451.84);
  draw p withcolor (0.00,0.00,0.00) dashed evenly scaled 1.82;
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 1.89;
  path p;
  p = (170.08, 607.75)
    --(510.24, 607.75);
  draw p withcolor (0.00,0.00,0.00) dashed evenly scaled 1.82;
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 1.89;
  path p;
  p = (184.25, 295.94)
    --(226.77, 678.61);
  draw p withcolor (0.00,0.00,0.00);
% Draw arrowhead type 0
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path arr;
  arr = (227.53, 668.85)
      --(226.52, 676.54)
      --(223.75, 669.29);
  draw arr withcolor (0.00,0.00,0.00);
% Draw arrowhead type 0
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path arr;
  arr = (183.43, 305.76)
      --(184.50, 298.02)
      --(187.21, 305.32);
  draw arr withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 1.89;
  path p;
  p = (170.08, 451.84)
    --(226.77, 678.61);
  draw p withcolor (0.00,0.00,0.00);
% Draw arrowhead type 0
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path arr;
  arr = (226.27, 668.85)
      --(226.27, 676.60)
      --(222.55, 669.73);
  draw arr withcolor (0.00,0.00,0.00);
% Draw arrowhead type 0
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path arr;
  arr = (170.58, 461.67)
      --(170.58, 453.86)
      --(174.24, 460.79);
  draw arr withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 1.89;
  path p;
  p = (184.25, 295.94)
    --(510.24, 451.84);
  draw p withcolor (0.00,0.00,0.00) dashed evenly scaled 1.82;
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 1.89;
  path p;
  p = (226.77, 338.46)
    --(372.28, 404.60);
  draw p withcolor (0.00,0.00,0.00);
% Draw arrowhead type 0
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path arr;
  arr = (364.28, 398.93)
      --(370.39, 403.72)
      --(362.71, 402.33);
  draw arr withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 1.89;
  path p;
  p = (245.67, 437.67)
    --(400.63, 437.67);
  draw p withcolor (0.00,0.00,0.00);
% Draw arrowhead type 0
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path arr;
  arr = (255.31, 439.56)
      --(247.75, 437.67)
      --(255.31, 435.78);
  draw arr withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 1.89;
  path p;
  p = (341.04, 354.65)
    --(486.55, 420.79);
  draw p withcolor (0.00,0.00,0.00);
% Draw arrowhead type 0
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path arr;
  arr = (348.98, 360.38)
      --(342.93, 355.53)
      --(350.55, 356.98);
  draw arr withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 1.89;
  path p;
  p = (170.08, 607.75)
    --(226.77, 678.61);
  draw p withcolor (0.00,0.00,0.00);
% Draw arrowhead type 0
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path arr;
  arr = (222.17, 669.92)
      --(225.45, 676.98)
      --(219.21, 672.31);
  draw arr withcolor (0.00,0.00,0.00);
% Draw arrowhead type 0
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path arr;
  arr = (174.61, 616.50)
      --(171.40, 609.39)
      --(177.57, 614.11);
  draw arr withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 1.89;
  path p;
  p = (333.67, 548.66)
    --(479.18, 482.52);
  draw p withcolor (0.00,0.00,0.00);
% Draw arrowhead type 0
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path arr;
  arr = (469.61, 484.85)
      --(477.29, 483.40)
      --(471.18, 488.25);
  draw arr withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 1.89;
  path p;
  p = (217.45, 573.61)
    --(362.96, 507.46);
  draw p withcolor (0.00,0.00,0.00);
% Draw arrowhead type 0
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path arr;
  arr = (226.96, 571.34)
      --(219.34, 572.72)
      --(225.39, 567.94);
  draw arr withcolor (0.00,0.00,0.00);
% End polyline object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt Sim\TeX etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((86.93,504.76))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt Sim\TeX etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((87.87,656.88))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt ${\rm C}_3$ etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((30.99,318.80))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt ${\rm C}_2$ etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((30.05,462.99))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt ${\rm C}_1$ etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((30.99,617.95))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt MatLab etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((93.54,359.24))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt uses etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((-14.17,324.28))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt uses etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((-14.17,466.02))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt uses etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((-14.17,621.92))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt contains etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((637.80,593.57))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt contains etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((637.80,451.84))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt contains etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((637.80,295.94))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt $(\overline{s}_N, \Delta_N)$ etex;
  picture p;
  p = q rotated 24.00;
  label.urt(p,((394.02,357.35))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt Truss/Stress $\mu$ etex;
  picture p;
  p = q rotated 24.00;
  label.urt(p,((222.05,346.02))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt Directory Service etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((196.98,720.19))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt $(\overline{s}_N, \Delta_N)$ etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((291.97,416.44))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt Truss/Deflection $\mu$ etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((298.14,318.61))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt Fin3D/Troot $\mu$ etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((295.75,629.48))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt $(\overline{s}_N, \Delta_N)$ etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((323.65,574.17))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt $(\overline{s}_N, \Delta_N)$ etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((331.15,261.92))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt Fin2D/Troot $\mu$ etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((296.69,474.52))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt $(\overline{s}_N, \Delta_N)$ etex;
  picture p;
  p = q rotated 336.00;
  label.urt(p,((239.94,510.80))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{18}{21.6}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt Fin3D/Troot $\mu$ etex;
  picture p;
  p = q rotated 336.00;
  label.urt(p,((354.77,501.42))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (515.72, 633.01)
    --(572.03, 633.01)
    --(572.03, 586.08)
    --(515.72, 586.08)
    --(515.72, 633.01)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.88,0.88,0.88) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (532.09, 586.08)
    --(555.59, 586.08)
    --(555.59, 583.75)
    --(532.09, 583.75)
    --(532.09, 586.08)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.50,0.50,0.50) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (529.76, 583.75)
    --(557.98, 583.75)
    --(557.98, 581.42)
    --(529.76, 581.42)
    --(529.76, 583.75)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.50,0.50,0.50) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (518.05, 564.98)
    --(569.70, 564.98)
    --(576.69, 553.20)
    --(576.69, 550.87)
    --(510.99, 550.87)
    --(510.99, 553.20)
    --(518.05, 564.98)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.88,0.88,0.88) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (516.85, 567.31)
    --(523.91, 567.31)
    --(523.91, 568.50)
    --(516.85, 568.50)
    --(516.85, 567.31)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (563.78, 567.31)
    --(570.83, 567.31)
    --(570.83, 568.50)
    --(563.78, 568.50)
    --(563.78, 567.31)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (516.85, 581.42)
    --(570.83, 581.42)
    --(570.83, 568.44)
    --(516.85, 568.44)
    --(516.85, 581.42)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.88,0.88,0.88) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (515.72, 335.37)
    --(572.03, 335.37)
    --(572.03, 288.44)
    --(515.72, 288.44)
    --(515.72, 335.37)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.88,0.88,0.88) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (532.09, 288.44)
    --(555.59, 288.44)
    --(555.59, 286.11)
    --(532.09, 286.11)
    --(532.09, 288.44)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.50,0.50,0.50) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (529.76, 286.11)
    --(557.98, 286.11)
    --(557.98, 283.78)
    --(529.76, 283.78)
    --(529.76, 286.11)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.50,0.50,0.50) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (518.05, 267.34)
    --(569.70, 267.34)
    --(576.69, 255.56)
    --(576.69, 253.23)
    --(510.99, 253.23)
    --(510.99, 255.56)
    --(518.05, 267.34)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.88,0.88,0.88) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (516.85, 269.67)
    --(523.91, 269.67)
    --(523.91, 270.87)
    --(516.85, 270.87)
    --(516.85, 269.67)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (563.78, 269.67)
    --(570.83, 269.67)
    --(570.83, 270.87)
    --(563.78, 270.87)
    --(563.78, 269.67)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (516.85, 283.78)
    --(570.83, 283.78)
    --(570.83, 270.80)
    --(516.85, 270.80)
    --(516.85, 283.78)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.88,0.88,0.88) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (66.65, 351.69)
    --(174.61, 351.69)
    --(174.61, 261.73)
    --(66.65, 261.73)
    --(66.65, 351.69)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.88,0.88,0.88) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (98.14, 261.73)
    --(143.12, 261.73)
    --(143.12, 257.20)
    --(98.14, 257.20)
    --(98.14, 261.73)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.50,0.50,0.50) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (93.61, 257.20)
    --(147.65, 257.20)
    --(147.65, 252.72)
    --(93.61, 252.72)
    --(93.61, 257.20)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.50,0.50,0.50) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (71.12, 252.72)
    --(170.14, 252.72)
    --(183.62, 230.24)
    --(183.62, 225.70)
    --(57.64, 225.70)
    --(57.64, 230.24)
    --(71.12, 252.72)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.88,0.88,0.88) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (516.66, 492.22)
    --(572.98, 492.22)
    --(572.98, 445.29)
    --(516.66, 445.29)
    --(516.66, 492.22)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.88,0.88,0.88) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (533.04, 445.29)
    --(556.54, 445.29)
    --(556.54, 442.96)
    --(533.04, 442.96)
    --(533.04, 445.29)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.50,0.50,0.50) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (530.71, 442.96)
    --(558.93, 442.96)
    --(558.93, 440.63)
    --(530.71, 440.63)
    --(530.71, 442.96)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.50,0.50,0.50) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (518.99, 424.19)
    --(570.65, 424.19)
    --(577.64, 412.41)
    --(577.64, 410.08)
    --(511.94, 410.08)
    --(511.94, 412.41)
    --(518.99, 424.19)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.88,0.88,0.88) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (517.80, 426.52)
    --(524.85, 426.52)
    --(524.85, 427.72)
    --(517.80, 427.72)
    --(517.80, 426.52)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (564.72, 426.52)
    --(571.78, 426.52)
    --(571.78, 427.72)
    --(564.72, 427.72)
    --(564.72, 426.52)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (517.80, 440.63)
    --(571.78, 440.63)
    --(571.78, 427.65)
    --(517.80, 427.65)
    --(517.80, 440.63)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.88,0.88,0.88) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=1;
  pickup pencircle scaled 0.00;
  path p;
  p = (560.31, 628.28)
    --(555.59, 628.28)
    --(555.59, 590.74)
    --(560.31, 590.74)
    --(560.31, 628.28)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.10white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=1;
  pickup pencircle scaled 0.00;
  path p;
  p = (555.59, 628.28)
    --(550.93, 628.28)
    --(550.93, 590.74)
    --(555.59, 590.74)
    --(555.59, 628.28)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.20white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=1;
  pickup pencircle scaled 0.00;
  path p;
  p = (550.93, 628.28)
    --(543.87, 628.28)
    --(543.87, 590.74)
    --(550.93, 590.74)
    --(550.93, 628.28)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.30white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=1;
  pickup pencircle scaled 0.00;
  path p;
  p = (543.87, 628.28)
    --(536.82, 628.28)
    --(536.82, 590.74)
    --(543.87, 590.74)
    --(543.87, 628.28)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.40white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=1;
  pickup pencircle scaled 0.00;
  path p;
  p = (536.82, 628.28)
    --(532.09, 628.28)
    --(532.09, 590.74)
    --(536.82, 590.74)
    --(536.82, 628.28)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.50white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=1;
  pickup pencircle scaled 0.00;
  path p;
  p = (532.09, 628.28)
    --(527.43, 628.28)
    --(527.43, 590.74)
    --(532.09, 590.74)
    --(532.09, 628.28)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.65white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=1;
  pickup pencircle scaled 0.00;
  path p;
  p = (527.43, 628.28)
    --(522.71, 628.28)
    --(522.71, 590.74)
    --(527.43, 590.74)
    --(527.43, 628.28)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.65white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=1;
  pickup pencircle scaled 0.00;
  path p;
  p = (564.98, 628.28)
    --(560.31, 628.28)
    --(560.31, 590.74)
    --(564.98, 590.74)
    --(564.98, 628.28)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (510.99, 553.20)
    --(510.99, 550.87)
    --(576.69, 550.87)
    --(576.69, 553.20);
  path f;
  f = p--cycle;
  fill f withcolor ((0.75,0.75,0.75) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (520.38, 562.58)
    --(515.72, 555.59)
    --(557.98, 555.59)
    --(555.59, 562.58)
    --(520.38, 562.58)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((1.00,1.00,1.00) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (567.37, 562.58)
    --(557.98, 562.58)
    --(560.31, 555.59)
    --(572.03, 555.59)
    --(567.37, 562.58)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((1.00,1.00,1.00) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (551.62, 576.57)
    --(551.62, 577.76)
    --(556.66, 577.76)
    --(556.66, 578.90)
    --(562.96, 578.90);
  draw p withcolor (0.50,0.50,0.50);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (556.66, 574.24)
    --(556.66, 576.57)
    --(562.96, 576.57);
  draw p withcolor (0.50,0.50,0.50);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (562.96, 577.76)
    --(568.06, 577.76);
  draw p withcolor (0.50,0.50,0.50);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (562.96, 578.90)
    --(562.96, 577.76);
  draw p withcolor (1.00,1.00,1.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (551.62, 576.57)
    --(556.66, 576.57);
  draw p withcolor (1.00,1.00,1.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (556.66, 574.24)
    --(562.96, 574.24)
    --(562.96, 576.57)
    --(568.06, 576.57)
    --(568.06, 577.76);
  draw p withcolor (1.00,1.00,1.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (564.22, 574.24)
    --(566.74, 574.24)
    --(566.74, 575.43);
  draw p withcolor (0.50,0.50,0.50);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (564.22, 574.24)
    --(564.22, 575.43)
    --(566.74, 575.43);
  draw p withcolor (1.00,1.00,1.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=1;
  pickup pencircle scaled 0.00;
  path p;
  p = (560.31, 330.65)
    --(555.59, 330.65)
    --(555.59, 293.10)
    --(560.31, 293.10)
    --(560.31, 330.65)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.10white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=1;
  pickup pencircle scaled 0.00;
  path p;
  p = (555.59, 330.65)
    --(550.93, 330.65)
    --(550.93, 293.10)
    --(555.59, 293.10)
    --(555.59, 330.65)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.20white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=1;
  pickup pencircle scaled 0.00;
  path p;
  p = (550.93, 330.65)
    --(543.87, 330.65)
    --(543.87, 293.10)
    --(550.93, 293.10)
    --(550.93, 330.65)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.30white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=1;
  pickup pencircle scaled 0.00;
  path p;
  p = (543.87, 330.65)
    --(536.82, 330.65)
    --(536.82, 293.10)
    --(543.87, 293.10)
    --(543.87, 330.65)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.40white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=1;
  pickup pencircle scaled 0.00;
  path p;
  p = (536.82, 330.65)
    --(532.09, 330.65)
    --(532.09, 293.10)
    --(536.82, 293.10)
    --(536.82, 330.65)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.50white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=1;
  pickup pencircle scaled 0.00;
  path p;
  p = (532.09, 330.65)
    --(527.43, 330.65)
    --(527.43, 293.10)
    --(532.09, 293.10)
    --(532.09, 330.65)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.65white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=1;
  pickup pencircle scaled 0.00;
  path p;
  p = (527.43, 330.65)
    --(522.71, 330.65)
    --(522.71, 293.10)
    --(527.43, 293.10)
    --(527.43, 330.65)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.65white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=1;
  pickup pencircle scaled 0.00;
  path p;
  p = (564.98, 330.65)
    --(560.31, 330.65)
    --(560.31, 293.10)
    --(564.98, 293.10)
    --(564.98, 330.65)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (510.99, 255.56)
    --(510.99, 253.23)
    --(576.69, 253.23)
    --(576.69, 255.56);
  path f;
  f = p--cycle;
  fill f withcolor ((0.75,0.75,0.75) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (520.38, 264.94)
    --(515.72, 257.95)
    --(557.98, 257.95)
    --(555.59, 264.94)
    --(520.38, 264.94)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((1.00,1.00,1.00) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (567.37, 264.94)
    --(557.98, 264.94)
    --(560.31, 257.95)
    --(572.03, 257.95)
    --(567.37, 264.94)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((1.00,1.00,1.00) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (551.62, 278.93)
    --(551.62, 280.13)
    --(556.66, 280.13)
    --(556.66, 281.26)
    --(562.96, 281.26);
  draw p withcolor (0.50,0.50,0.50);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (556.66, 276.60)
    --(556.66, 278.93)
    --(562.96, 278.93);
  draw p withcolor (0.50,0.50,0.50);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (562.96, 280.13)
    --(568.06, 280.13);
  draw p withcolor (0.50,0.50,0.50);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (562.96, 281.26)
    --(562.96, 280.13);
  draw p withcolor (1.00,1.00,1.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (551.62, 278.93)
    --(556.66, 278.93);
  draw p withcolor (1.00,1.00,1.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (556.66, 276.60)
    --(562.96, 276.60)
    --(562.96, 278.93)
    --(568.06, 278.93)
    --(568.06, 280.13);
  draw p withcolor (1.00,1.00,1.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (564.22, 276.60)
    --(566.74, 276.60)
    --(566.74, 277.80);
  draw p withcolor (0.50,0.50,0.50);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (564.22, 276.60)
    --(564.22, 277.80)
    --(566.74, 277.80);
  draw p withcolor (1.00,1.00,1.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (57.64, 230.24)
    --(57.64, 225.70)
    --(183.62, 225.70)
    --(183.62, 230.24);
  path f;
  f = p--cycle;
  fill f withcolor ((0.75,0.75,0.75) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (75.65, 248.19)
    --(66.65, 234.71)
    --(147.65, 234.71)
    --(143.12, 248.19)
    --(75.65, 248.19)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((1.00,1.00,1.00) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (165.61, 248.19)
    --(147.65, 248.19)
    --(152.13, 234.71)
    --(174.61, 234.71)
    --(165.61, 248.19)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((1.00,1.00,1.00) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.00;
  path p;
  p = (80.13, 342.74)
    --(161.13, 342.74)
    --(161.13, 270.74)
    --(80.13, 270.74)
    --(80.13, 342.74)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=1;
  pickup pencircle scaled 0.00;
  path p;
  p = (561.26, 487.50)
    --(556.54, 487.50)
    --(556.54, 449.95)
    --(561.26, 449.95)
    --(561.26, 487.50)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.10white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=1;
  pickup pencircle scaled 0.00;
  path p;
  p = (556.54, 487.50)
    --(551.87, 487.50)
    --(551.87, 449.95)
    --(556.54, 449.95)
    --(556.54, 487.50)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.20white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=1;
  pickup pencircle scaled 0.00;
  path p;
  p = (551.87, 487.50)
    --(544.82, 487.50)
    --(544.82, 449.95)
    --(551.87, 449.95)
    --(551.87, 487.50)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.30white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=1;
  pickup pencircle scaled 0.00;
  path p;
  p = (544.82, 487.50)
    --(537.76, 487.50)
    --(537.76, 449.95)
    --(544.82, 449.95)
    --(544.82, 487.50)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.40white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=1;
  pickup pencircle scaled 0.00;
  path p;
  p = (537.76, 487.50)
    --(533.04, 487.50)
    --(533.04, 449.95)
    --(537.76, 449.95)
    --(537.76, 487.50)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.50white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=1;
  pickup pencircle scaled 0.00;
  path p;
  p = (533.04, 487.50)
    --(528.38, 487.50)
    --(528.38, 449.95)
    --(533.04, 449.95)
    --(533.04, 487.50)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.65white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=1;
  pickup pencircle scaled 0.00;
  path p;
  p = (528.38, 487.50)
    --(523.65, 487.50)
    --(523.65, 449.95)
    --(528.38, 449.95)
    --(528.38, 487.50)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.65white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=1;
  pickup pencircle scaled 0.00;
  path p;
  p = (565.92, 487.50)
    --(561.26, 487.50)
    --(561.26, 449.95)
    --(565.92, 449.95)
    --(565.92, 487.50)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.25,0.25) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (511.94, 412.41)
    --(511.94, 410.08)
    --(577.64, 410.08)
    --(577.64, 412.41);
  path f;
  f = p--cycle;
  fill f withcolor ((0.75,0.75,0.75) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (521.32, 421.80)
    --(516.66, 414.80)
    --(558.93, 414.80)
    --(556.54, 421.80)
    --(521.32, 421.80)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((1.00,1.00,1.00) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (568.31, 421.80)
    --(558.93, 421.80)
    --(561.26, 414.80)
    --(572.98, 414.80)
    --(568.31, 421.80)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((1.00,1.00,1.00) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (552.57, 435.78)
    --(552.57, 436.98)
    --(557.61, 436.98)
    --(557.61, 438.11)
    --(563.91, 438.11);
  draw p withcolor (0.50,0.50,0.50);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (557.61, 433.45)
    --(557.61, 435.78)
    --(563.91, 435.78);
  draw p withcolor (0.50,0.50,0.50);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (563.91, 436.98)
    --(569.01, 436.98);
  draw p withcolor (0.50,0.50,0.50);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (563.91, 438.11)
    --(563.91, 436.98);
  draw p withcolor (1.00,1.00,1.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (552.57, 435.78)
    --(557.61, 435.78);
  draw p withcolor (1.00,1.00,1.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (557.61, 433.45)
    --(563.91, 433.45)
    --(563.91, 435.78)
    --(569.01, 435.78)
    --(569.01, 436.98);
  draw p withcolor (1.00,1.00,1.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (565.17, 433.45)
    --(567.69, 433.45)
    --(567.69, 434.65);
  draw p withcolor (0.50,0.50,0.50);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (565.17, 433.45)
    --(565.17, 434.65)
    --(567.69, 434.65);
  draw p withcolor (1.00,1.00,1.00);
% End polyline object
% Begin ellipse object
  pickup pencircle scaled 0.94;
  path c;
  c = fullcircle scaled 1.26 yscaled 1.00
         rotated 0.00 shifted (243.84,690.83);
  draw c withcolor (0.00,0.00,0.00);
% End ellipse object
% Begin ellipse object
  pickup pencircle scaled 0.94;
  path c;
  c = fullcircle scaled 1.26 yscaled 1.00
         rotated 0.00 shifted (253.04,690.83);
  draw c withcolor (0.00,0.00,0.00);
% End ellipse object
% Begin polyline object
  linecap:=1;
  linejoin:=2;
  pickup pencircle scaled 0.94;
  path p;
  p = (244.72, 691.91)
    --(244.79, 691.59)
    --(244.85, 691.15)
    --(244.91, 690.58)
    --(245.04, 690.14)
    --(245.17, 689.76)
    --(245.29, 689.57)
    --(245.42, 689.39)
    --(245.61, 689.32)
    --(245.73, 689.26)
    --(245.80, 689.26)
    --(245.86, 689.26)
    --(246.30, 689.26)
    --(247.18, 689.26)
    --(248.44, 689.26)
    --(249.70, 689.26)
    --(250.58, 689.26)
    --(251.02, 689.26)
    --(251.09, 689.26)
    --(251.15, 689.26)
    --(251.28, 689.32)
    --(251.46, 689.39)
    --(251.59, 689.57)
    --(251.72, 689.76)
    --(251.84, 690.14)
    --(251.97, 690.58)
    --(252.09, 691.15)
    --(252.16, 691.59)
    --(252.16, 691.91)
    --(252.16, 692.09)
    --(252.16, 692.16)
    --(252.09, 692.28)
    --(252.03, 692.35)
    --(251.97, 692.41)
    --(251.91, 692.41)
    --(251.78, 692.47)
    --(251.72, 692.47)
    --(251.65, 692.47)
    --(251.40, 692.47)
    --(250.71, 692.47)
    --(249.70, 692.47)
    --(248.44, 692.47)
    --(247.18, 692.47)
    --(246.17, 692.47)
    --(245.48, 692.47)
    --(245.23, 692.47)
    --(245.17, 692.47)
    --(245.10, 692.47)
    --(244.98, 692.41)
    --(244.91, 692.41)
    --(244.85, 692.35)
    --(244.79, 692.28)
    --(244.79, 692.16)
    --(244.72, 692.09)
    --(244.72, 691.91)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.53,0.51,0.53) + 0.00white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (520.38, 590.74)
    --(520.44, 590.80)
    --(520.82, 591.18)
    --(521.57, 591.94)
    --(522.27, 592.69)
    --(522.65, 593.07)
    --(522.71, 593.13)
    --(522.71, 593.32)
    --(522.71, 593.70)
    --(522.71, 594.39)
    --(522.71, 595.46)
    --(522.71, 596.98)
    --(522.71, 598.87)
    --(522.71, 601.13)
    --(522.71, 603.72)
    --(522.71, 606.55)
    --(522.71, 609.57)
    --(522.71, 612.54)
    --(522.71, 615.37)
    --(522.71, 617.95)
    --(522.71, 620.22)
    --(522.71, 622.11)
    --(522.71, 623.62)
    --(522.71, 624.69)
    --(522.71, 625.39)
    --(522.71, 625.76)
    --(522.71, 625.95)
    --(522.65, 626.02)
    --(522.27, 626.39)
    --(521.51, 627.15)
    --(520.82, 627.84)
    --(520.44, 628.22)
    --(520.38, 628.28)
    --(520.38, 628.16)
    --(520.38, 627.78)
    --(520.38, 627.15)
    --(520.38, 626.14)
    --(520.38, 624.82)
    --(520.38, 623.06)
    --(520.38, 620.91)
    --(520.38, 618.39)
    --(520.38, 615.62)
    --(520.38, 612.60)
    --(520.38, 609.51)
    --(520.38, 606.43)
    --(520.38, 603.40)
    --(520.38, 600.63)
    --(520.38, 598.11)
    --(520.38, 595.97)
    --(520.38, 594.20)
    --(520.38, 592.88)
    --(520.38, 591.87)
    --(520.38, 591.24)
    --(520.38, 590.87)
    --(520.38, 590.74)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.75,0.75,0.75) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (522.71, 593.13)
    --(522.83, 593.13)
    --(523.02, 593.13)
    --(523.46, 593.07)
    --(524.09, 593.01)
    --(524.85, 592.94)
    --(525.86, 592.88)
    --(526.99, 592.82)
    --(528.25, 592.69)
    --(529.57, 592.63)
    --(531.02, 592.57)
    --(532.54, 592.44)
    --(534.11, 592.38)
    --(535.87, 592.31)
    --(537.76, 592.25)
    --(539.84, 592.25)
    --(542.05, 592.19)
    --(544.38, 592.19)
    --(546.71, 592.19)
    --(548.85, 592.25)
    --(550.87, 592.25)
    --(552.69, 592.31)
    --(554.39, 592.38)
    --(555.91, 592.44)
    --(557.35, 592.57)
    --(558.68, 592.63)
    --(559.94, 592.69)
    --(561.07, 592.82)
    --(562.08, 592.88)
    --(563.02, 592.94)
    --(563.72, 593.01)
    --(564.28, 593.07)
    --(564.66, 593.13)
    --(564.85, 593.13)
    --(564.98, 593.13)
    --(565.04, 593.07)
    --(565.42, 592.69)
    --(566.17, 591.94)
    --(566.93, 591.18)
    --(567.31, 590.80)
    --(567.37, 590.74)
    --(567.24, 590.74)
    --(566.87, 590.74)
    --(566.24, 590.74)
    --(565.23, 590.74)
    --(563.84, 590.74)
    --(562.08, 590.74)
    --(559.87, 590.74)
    --(557.23, 590.74)
    --(554.27, 590.74)
    --(550.93, 590.74)
    --(547.46, 590.74)
    --(543.87, 590.74)
    --(540.28, 590.74)
    --(536.82, 590.74)
    --(533.48, 590.74)
    --(530.52, 590.74)
    --(527.87, 590.74)
    --(525.67, 590.74)
    --(523.91, 590.74)
    --(522.52, 590.74)
    --(521.51, 590.74)
    --(520.88, 590.74)
    --(520.50, 590.74)
    --(520.38, 590.74)
    --(520.44, 590.80)
    --(520.82, 591.18)
    --(521.57, 591.94)
    --(522.27, 592.69)
    --(522.65, 593.07)
    --(522.71, 593.13)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((1.00,1.00,1.00) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (567.37, 590.74)
    --(567.31, 590.80)
    --(566.93, 591.18)
    --(566.17, 591.94)
    --(565.42, 592.69)
    --(565.04, 593.07)
    --(564.98, 593.13)
    --(564.98, 593.32)
    --(564.98, 593.70)
    --(564.98, 594.39)
    --(564.98, 595.46)
    --(564.98, 596.98)
    --(564.98, 598.87)
    --(564.98, 601.13)
    --(564.98, 603.72)
    --(564.98, 606.55)
    --(564.98, 609.57)
    --(564.98, 612.54)
    --(564.98, 615.37)
    --(564.98, 617.95)
    --(564.98, 620.22)
    --(564.98, 622.11)
    --(564.98, 623.62)
    --(564.98, 624.69)
    --(564.98, 625.39)
    --(564.98, 625.76)
    --(564.98, 625.95)
    --(565.04, 626.02)
    --(565.42, 626.39)
    --(566.17, 627.15)
    --(566.93, 627.84)
    --(567.31, 628.22)
    --(567.37, 628.28)
    --(567.37, 628.16)
    --(567.37, 627.78)
    --(567.37, 627.15)
    --(567.37, 626.14)
    --(567.37, 624.82)
    --(567.37, 623.06)
    --(567.37, 620.91)
    --(567.37, 618.39)
    --(567.37, 615.62)
    --(567.37, 612.60)
    --(567.37, 609.51)
    --(567.37, 606.43)
    --(567.37, 603.40)
    --(567.37, 600.63)
    --(567.37, 598.11)
    --(567.37, 595.97)
    --(567.37, 594.20)
    --(567.37, 592.88)
    --(567.37, 591.87)
    --(567.37, 591.24)
    --(567.37, 590.87)
    --(567.37, 590.74)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.75,0.75,0.75) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (522.71, 625.95)
    --(522.83, 625.95)
    --(523.02, 625.95)
    --(523.46, 626.02)
    --(524.09, 626.08)
    --(524.85, 626.14)
    --(525.86, 626.20)
    --(526.99, 626.33)
    --(528.25, 626.39)
    --(529.57, 626.46)
    --(531.02, 626.58)
    --(532.54, 626.65)
    --(534.11, 626.71)
    --(535.87, 626.83)
    --(537.76, 626.90)
    --(539.84, 626.90)
    --(542.05, 626.96)
    --(544.38, 626.96)
    --(546.71, 626.96)
    --(548.85, 626.90)
    --(550.87, 626.90)
    --(552.69, 626.83)
    --(554.39, 626.71)
    --(555.91, 626.65)
    --(557.35, 626.58)
    --(558.68, 626.46)
    --(559.94, 626.39)
    --(561.07, 626.33)
    --(562.08, 626.20)
    --(563.02, 626.14)
    --(563.72, 626.08)
    --(564.28, 626.02)
    --(564.66, 625.95)
    --(564.85, 625.95)
    --(564.98, 625.95)
    --(565.04, 626.02)
    --(565.42, 626.39)
    --(566.17, 627.15)
    --(566.93, 627.84)
    --(567.31, 628.22)
    --(567.37, 628.28)
    --(567.24, 628.28)
    --(566.87, 628.28)
    --(566.24, 628.28)
    --(565.23, 628.28)
    --(563.84, 628.28)
    --(562.08, 628.28)
    --(559.87, 628.28)
    --(557.23, 628.28)
    --(554.27, 628.28)
    --(550.93, 628.28)
    --(547.46, 628.28)
    --(543.87, 628.28)
    --(540.28, 628.28)
    --(536.82, 628.28)
    --(533.48, 628.28)
    --(530.52, 628.28)
    --(527.87, 628.28)
    --(525.67, 628.28)
    --(523.91, 628.28)
    --(522.52, 628.28)
    --(521.51, 628.28)
    --(520.88, 628.28)
    --(520.50, 628.28)
    --(520.38, 628.28)
    --(520.44, 628.22)
    --(520.82, 627.84)
    --(521.57, 627.09)
    --(522.27, 626.39)
    --(522.65, 626.02)
    --(522.71, 625.95)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.50,0.50,0.50) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (515.72, 555.59)
    --(557.98, 555.59)
    --(557.98, 554.71)
    --(515.72, 554.71)
    --(515.72, 555.59)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.50,0.50,0.50) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (519.50, 562.58)
    --(515.72, 555.59)
    --(522.71, 555.59)
    --(523.72, 557.92)
    --(522.08, 557.92)
    --(523.02, 560.25)
    --(521.51, 560.25)
    --(522.71, 562.58)
    --(519.50, 562.58)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.75,0.75,0.75) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (557.98, 555.59)
    --(550.93, 555.59)
    --(550.93, 557.92)
    --(552.44, 557.92)
    --(552.19, 560.25)
    --(553.26, 560.25)
    --(553.57, 560.25)
    --(553.13, 562.58)
    --(556.03, 562.58)
    --(557.98, 555.59)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.75,0.75,0.75) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (572.03, 555.59)
    --(569.89, 555.59)
    --(567.37, 561.20)
    --(558.36, 561.20)
    --(557.92, 562.58)
    --(568.31, 562.58)
    --(572.03, 555.59)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.75,0.75,0.75) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (560.31, 555.59)
    --(572.03, 555.59)
    --(572.03, 554.77)
    --(560.31, 554.77)
    --(560.31, 555.59)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.50,0.50,0.50) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (560.31, 555.53)
    --(560.25, 554.71)
    --(557.92, 561.89)
    --(557.92, 562.58)
    --(560.31, 555.53)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.50,0.50,0.50) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (520.38, 293.10)
    --(520.44, 293.17)
    --(520.82, 293.54)
    --(521.57, 294.30)
    --(522.27, 295.06)
    --(522.65, 295.43)
    --(522.71, 295.50)
    --(522.71, 295.69)
    --(522.71, 296.06)
    --(522.71, 296.76)
    --(522.71, 297.83)
    --(522.71, 299.34)
    --(522.71, 301.23)
    --(522.71, 303.50)
    --(522.71, 306.08)
    --(522.71, 308.91)
    --(522.71, 311.94)
    --(522.71, 314.90)
    --(522.71, 317.73)
    --(522.71, 320.31)
    --(522.71, 322.58)
    --(522.71, 324.47)
    --(522.71, 325.98)
    --(522.71, 327.06)
    --(522.71, 327.75)
    --(522.71, 328.13)
    --(522.71, 328.31)
    --(522.65, 328.38)
    --(522.27, 328.76)
    --(521.51, 329.51)
    --(520.82, 330.20)
    --(520.44, 330.58)
    --(520.38, 330.65)
    --(520.38, 330.52)
    --(520.38, 330.14)
    --(520.38, 329.51)
    --(520.38, 328.50)
    --(520.38, 327.18)
    --(520.38, 325.42)
    --(520.38, 323.28)
    --(520.38, 320.76)
    --(520.38, 317.98)
    --(520.38, 314.96)
    --(520.38, 311.87)
    --(520.38, 308.79)
    --(520.38, 305.76)
    --(520.38, 302.99)
    --(520.38, 300.47)
    --(520.38, 298.33)
    --(520.38, 296.57)
    --(520.38, 295.24)
    --(520.38, 294.24)
    --(520.38, 293.61)
    --(520.38, 293.23)
    --(520.38, 293.10)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.75,0.75,0.75) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (522.71, 295.50)
    --(522.83, 295.50)
    --(523.02, 295.50)
    --(523.46, 295.43)
    --(524.09, 295.37)
    --(524.85, 295.31)
    --(525.86, 295.24)
    --(526.99, 295.18)
    --(528.25, 295.06)
    --(529.57, 294.99)
    --(531.02, 294.93)
    --(532.54, 294.80)
    --(534.11, 294.74)
    --(535.87, 294.68)
    --(537.76, 294.61)
    --(539.84, 294.61)
    --(542.05, 294.55)
    --(544.38, 294.55)
    --(546.71, 294.55)
    --(548.85, 294.61)
    --(550.87, 294.61)
    --(552.69, 294.68)
    --(554.39, 294.74)
    --(555.91, 294.80)
    --(557.35, 294.93)
    --(558.68, 294.99)
    --(559.94, 295.06)
    --(561.07, 295.18)
    --(562.08, 295.24)
    --(563.02, 295.31)
    --(563.72, 295.37)
    --(564.28, 295.43)
    --(564.66, 295.50)
    --(564.85, 295.50)
    --(564.98, 295.50)
    --(565.04, 295.43)
    --(565.42, 295.06)
    --(566.17, 294.30)
    --(566.93, 293.54)
    --(567.31, 293.17)
    --(567.37, 293.10)
    --(567.24, 293.10)
    --(566.87, 293.10)
    --(566.24, 293.10)
    --(565.23, 293.10)
    --(563.84, 293.10)
    --(562.08, 293.10)
    --(559.87, 293.10)
    --(557.23, 293.10)
    --(554.27, 293.10)
    --(550.93, 293.10)
    --(547.46, 293.10)
    --(543.87, 293.10)
    --(540.28, 293.10)
    --(536.82, 293.10)
    --(533.48, 293.10)
    --(530.52, 293.10)
    --(527.87, 293.10)
    --(525.67, 293.10)
    --(523.91, 293.10)
    --(522.52, 293.10)
    --(521.51, 293.10)
    --(520.88, 293.10)
    --(520.50, 293.10)
    --(520.38, 293.10)
    --(520.44, 293.17)
    --(520.82, 293.54)
    --(521.57, 294.30)
    --(522.27, 295.06)
    --(522.65, 295.43)
    --(522.71, 295.50)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((1.00,1.00,1.00) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (567.37, 293.10)
    --(567.31, 293.17)
    --(566.93, 293.54)
    --(566.17, 294.30)
    --(565.42, 295.06)
    --(565.04, 295.43)
    --(564.98, 295.50)
    --(564.98, 295.69)
    --(564.98, 296.06)
    --(564.98, 296.76)
    --(564.98, 297.83)
    --(564.98, 299.34)
    --(564.98, 301.23)
    --(564.98, 303.50)
    --(564.98, 306.08)
    --(564.98, 308.91)
    --(564.98, 311.94)
    --(564.98, 314.90)
    --(564.98, 317.73)
    --(564.98, 320.31)
    --(564.98, 322.58)
    --(564.98, 324.47)
    --(564.98, 325.98)
    --(564.98, 327.06)
    --(564.98, 327.75)
    --(564.98, 328.13)
    --(564.98, 328.31)
    --(565.04, 328.38)
    --(565.42, 328.76)
    --(566.17, 329.51)
    --(566.93, 330.20)
    --(567.31, 330.58)
    --(567.37, 330.65)
    --(567.37, 330.52)
    --(567.37, 330.14)
    --(567.37, 329.51)
    --(567.37, 328.50)
    --(567.37, 327.18)
    --(567.37, 325.42)
    --(567.37, 323.28)
    --(567.37, 320.76)
    --(567.37, 317.98)
    --(567.37, 314.96)
    --(567.37, 311.87)
    --(567.37, 308.79)
    --(567.37, 305.76)
    --(567.37, 302.99)
    --(567.37, 300.47)
    --(567.37, 298.33)
    --(567.37, 296.57)
    --(567.37, 295.24)
    --(567.37, 294.24)
    --(567.37, 293.61)
    --(567.37, 293.23)
    --(567.37, 293.10)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.75,0.75,0.75) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (522.71, 328.31)
    --(522.83, 328.31)
    --(523.02, 328.31)
    --(523.46, 328.38)
    --(524.09, 328.44)
    --(524.85, 328.50)
    --(525.86, 328.57)
    --(526.99, 328.69)
    --(528.25, 328.76)
    --(529.57, 328.82)
    --(531.02, 328.94)
    --(532.54, 329.01)
    --(534.11, 329.07)
    --(535.87, 329.20)
    --(537.76, 329.26)
    --(539.84, 329.26)
    --(542.05, 329.32)
    --(544.38, 329.32)
    --(546.71, 329.32)
    --(548.85, 329.26)
    --(550.87, 329.26)
    --(552.69, 329.20)
    --(554.39, 329.07)
    --(555.91, 329.01)
    --(557.35, 328.94)
    --(558.68, 328.82)
    --(559.94, 328.76)
    --(561.07, 328.69)
    --(562.08, 328.57)
    --(563.02, 328.50)
    --(563.72, 328.44)
    --(564.28, 328.38)
    --(564.66, 328.31)
    --(564.85, 328.31)
    --(564.98, 328.31)
    --(565.04, 328.38)
    --(565.42, 328.76)
    --(566.17, 329.51)
    --(566.93, 330.20)
    --(567.31, 330.58)
    --(567.37, 330.65)
    --(567.24, 330.65)
    --(566.87, 330.65)
    --(566.24, 330.65)
    --(565.23, 330.65)
    --(563.84, 330.65)
    --(562.08, 330.65)
    --(559.87, 330.65)
    --(557.23, 330.65)
    --(554.27, 330.65)
    --(550.93, 330.65)
    --(547.46, 330.65)
    --(543.87, 330.65)
    --(540.28, 330.65)
    --(536.82, 330.65)
    --(533.48, 330.65)
    --(530.52, 330.65)
    --(527.87, 330.65)
    --(525.67, 330.65)
    --(523.91, 330.65)
    --(522.52, 330.65)
    --(521.51, 330.65)
    --(520.88, 330.65)
    --(520.50, 330.65)
    --(520.38, 330.65)
    --(520.44, 330.58)
    --(520.82, 330.20)
    --(521.57, 329.45)
    --(522.27, 328.76)
    --(522.65, 328.38)
    --(522.71, 328.31)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.50,0.50,0.50) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (515.72, 257.95)
    --(557.98, 257.95)
    --(557.98, 257.07)
    --(515.72, 257.07)
    --(515.72, 257.95)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.50,0.50,0.50) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (519.50, 264.94)
    --(515.72, 257.95)
    --(522.71, 257.95)
    --(523.72, 260.28)
    --(522.08, 260.28)
    --(523.02, 262.61)
    --(521.51, 262.61)
    --(522.71, 264.94)
    --(519.50, 264.94)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.75,0.75,0.75) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (557.98, 257.95)
    --(550.93, 257.95)
    --(550.93, 260.28)
    --(552.44, 260.28)
    --(552.19, 262.61)
    --(553.26, 262.61)
    --(553.57, 262.61)
    --(553.13, 264.94)
    --(556.03, 264.94)
    --(557.98, 257.95)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.75,0.75,0.75) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (572.03, 257.95)
    --(569.89, 257.95)
    --(567.37, 263.56)
    --(558.36, 263.56)
    --(557.92, 264.94)
    --(568.31, 264.94)
    --(572.03, 257.95)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.75,0.75,0.75) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (560.31, 257.95)
    --(572.03, 257.95)
    --(572.03, 257.13)
    --(560.31, 257.13)
    --(560.31, 257.95)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.50,0.50,0.50) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (560.31, 257.89)
    --(560.25, 257.07)
    --(557.92, 264.25)
    --(557.92, 264.94)
    --(560.31, 257.89)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.50,0.50,0.50) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (66.65, 234.71)
    --(147.65, 234.71)
    --(147.65, 233.07)
    --(66.65, 233.07)
    --(66.65, 234.71)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.50,0.50,0.50) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (73.95, 248.19)
    --(66.65, 234.71)
    --(80.13, 234.71)
    --(82.02, 239.24)
    --(78.87, 239.24)
    --(80.69, 243.72)
    --(77.80, 243.72)
    --(80.13, 248.19)
    --(73.95, 248.19)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.75,0.75,0.75) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (147.65, 234.71)
    --(134.11, 234.71)
    --(134.11, 239.24)
    --(137.07, 239.24)
    --(136.57, 243.72)
    --(138.65, 243.72)
    --(139.28, 243.72)
    --(138.39, 248.19)
    --(144.00, 248.19)
    --(147.65, 234.71)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.75,0.75,0.75) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (174.61, 234.71)
    --(170.52, 234.71)
    --(165.61, 245.54)
    --(148.47, 245.54)
    --(147.59, 248.19)
    --(167.50, 248.19)
    --(174.61, 234.71)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.75,0.75,0.75) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (152.13, 234.71)
    --(174.61, 234.71)
    --(174.61, 233.13)
    --(152.13, 233.13)
    --(152.13, 234.71)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.50,0.50,0.50) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (152.13, 234.65)
    --(152.06, 233.01)
    --(147.53, 246.87)
    --(147.59, 248.13)
    --(152.13, 234.65)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.50,0.50,0.50) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (80.13, 275.21)
    --(80.19, 275.21)
    --(80.44, 275.21)
    --(80.76, 275.15)
    --(81.32, 275.09)
    --(82.08, 275.02)
    --(83.09, 274.96)
    --(84.22, 274.90)
    --(85.61, 274.77)
    --(87.18, 274.65)
    --(88.82, 274.52)
    --(90.65, 274.39)
    --(92.54, 274.27)
    --(94.55, 274.14)
    --(96.63, 274.02)
    --(98.83, 273.95)
    --(101.10, 273.83)
    --(103.56, 273.70)
    --(106.14, 273.64)
    --(108.85, 273.57)
    --(111.81, 273.51)
    --(114.90, 273.45)
    --(118.17, 273.39)
    --(121.51, 273.39)
    --(124.85, 273.39)
    --(128.06, 273.45)
    --(131.15, 273.51)
    --(133.98, 273.57)
    --(136.69, 273.64)
    --(139.15, 273.70)
    --(141.48, 273.83)
    --(143.69, 273.95)
    --(145.76, 274.02)
    --(147.72, 274.14)
    --(149.61, 274.27)
    --(151.37, 274.39)
    --(153.07, 274.52)
    --(154.65, 274.65)
    --(156.03, 274.77)
    --(157.35, 274.90)
    --(158.43, 274.96)
    --(159.31, 275.02)
    --(160.00, 275.09)
    --(160.57, 275.15)
    --(160.88, 275.21)
    --(161.07, 275.21)
    --(161.13, 275.21)
    --(161.20, 275.15)
    --(161.39, 274.96)
    --(161.83, 274.52)
    --(162.52, 273.83)
    --(163.40, 272.94)
    --(164.22, 272.13)
    --(164.91, 271.43)
    --(165.35, 270.99)
    --(165.54, 270.80)
    --(165.61, 270.74)
    --(165.48, 270.74)
    --(165.23, 270.74)
    --(164.72, 270.74)
    --(163.97, 270.74)
    --(162.90, 270.74)
    --(161.45, 270.74)
    --(159.56, 270.74)
    --(157.29, 270.74)
    --(154.58, 270.74)
    --(151.50, 270.74)
    --(147.97, 270.74)
    --(144.06, 270.74)
    --(139.78, 270.74)
    --(135.24, 270.74)
    --(130.52, 270.74)
    --(125.61, 270.74)
    --(120.63, 270.74)
    --(115.65, 270.74)
    --(110.74, 270.74)
    --(106.02, 270.74)
    --(101.48, 270.74)
    --(97.20, 270.74)
    --(93.29, 270.74)
    --(89.76, 270.74)
    --(86.68, 270.74)
    --(83.97, 270.74)
    --(81.70, 270.74)
    --(79.81, 270.74)
    --(78.36, 270.74)
    --(77.29, 270.74)
    --(76.54, 270.74)
    --(76.03, 270.74)
    --(75.78, 270.74)
    --(75.65, 270.74)
    --(75.72, 270.80)
    --(75.91, 270.99)
    --(76.35, 271.43)
    --(77.04, 272.13)
    --(77.92, 273.01)
    --(78.74, 273.83)
    --(79.43, 274.52)
    --(79.87, 274.96)
    --(80.06, 275.15)
    --(80.13, 275.21)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((1.00,1.00,1.00) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (165.61, 270.74)
    --(165.54, 270.80)
    --(165.35, 270.99)
    --(164.91, 271.43)
    --(164.22, 272.13)
    --(163.34, 273.01)
    --(162.52, 273.83)
    --(161.83, 274.52)
    --(161.39, 274.96)
    --(161.20, 275.15)
    --(161.13, 275.21)
    --(161.13, 275.34)
    --(161.13, 275.65)
    --(161.13, 276.22)
    --(161.13, 277.10)
    --(161.13, 278.36)
    --(161.13, 280.06)
    --(161.13, 282.14)
    --(161.13, 284.66)
    --(161.13, 287.56)
    --(161.13, 290.90)
    --(161.13, 294.55)
    --(161.13, 298.46)
    --(161.13, 302.55)
    --(161.13, 306.71)
    --(161.13, 310.87)
    --(161.13, 314.96)
    --(161.13, 318.87)
    --(161.13, 322.52)
    --(161.13, 325.86)
    --(161.13, 328.76)
    --(161.13, 331.28)
    --(161.13, 333.35)
    --(161.13, 335.06)
    --(161.13, 336.31)
    --(161.13, 337.20)
    --(161.13, 337.76)
    --(161.13, 338.08)
    --(161.13, 338.20)
    --(161.20, 338.27)
    --(161.39, 338.46)
    --(161.83, 338.90)
    --(162.52, 339.59)
    --(163.40, 340.47)
    --(164.22, 341.35)
    --(164.91, 342.05)
    --(165.35, 342.49)
    --(165.54, 342.68)
    --(165.61, 342.74)
    --(165.61, 342.61)
    --(165.61, 342.30)
    --(165.61, 341.80)
    --(165.61, 340.91)
    --(165.61, 339.72)
    --(165.61, 338.08)
    --(165.61, 336.06)
    --(165.61, 333.61)
    --(165.61, 330.77)
    --(165.61, 327.50)
    --(165.61, 323.84)
    --(165.61, 319.81)
    --(165.61, 315.59)
    --(165.61, 311.24)
    --(165.61, 306.71)
    --(165.61, 302.24)
    --(165.61, 297.89)
    --(165.61, 293.67)
    --(165.61, 289.64)
    --(165.61, 285.98)
    --(165.61, 282.71)
    --(165.61, 279.87)
    --(165.61, 277.42)
    --(165.61, 275.40)
    --(165.61, 273.76)
    --(165.61, 272.57)
    --(165.61, 271.69)
    --(165.61, 271.18)
    --(165.61, 270.87)
    --(165.61, 270.74)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.75,0.75,0.75) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (75.65, 270.74)
    --(75.72, 270.80)
    --(75.91, 270.99)
    --(76.35, 271.43)
    --(77.04, 272.13)
    --(77.92, 273.01)
    --(78.74, 273.83)
    --(79.43, 274.52)
    --(79.87, 274.96)
    --(80.06, 275.15)
    --(80.13, 275.21)
    --(80.13, 275.34)
    --(80.13, 275.65)
    --(80.13, 276.22)
    --(80.13, 277.10)
    --(80.13, 278.36)
    --(80.13, 280.06)
    --(80.13, 282.14)
    --(80.13, 284.66)
    --(80.13, 287.56)
    --(80.13, 290.90)
    --(80.13, 294.55)
    --(80.13, 298.46)
    --(80.13, 302.55)
    --(80.13, 306.71)
    --(80.13, 310.87)
    --(80.13, 314.96)
    --(80.13, 318.87)
    --(80.13, 322.52)
    --(80.13, 325.86)
    --(80.13, 328.76)
    --(80.13, 331.28)
    --(80.13, 333.35)
    --(80.13, 335.06)
    --(80.13, 336.31)
    --(80.13, 337.20)
    --(80.13, 337.76)
    --(80.13, 338.08)
    --(80.13, 338.20)
    --(80.06, 338.27)
    --(79.87, 338.46)
    --(79.43, 338.90)
    --(78.74, 339.59)
    --(77.86, 340.47)
    --(77.04, 341.35)
    --(76.35, 342.05)
    --(75.91, 342.49)
    --(75.72, 342.68)
    --(75.65, 342.74)
    --(75.65, 342.61)
    --(75.65, 342.30)
    --(75.65, 341.80)
    --(75.65, 340.91)
    --(75.65, 339.72)
    --(75.65, 338.08)
    --(75.65, 336.06)
    --(75.65, 333.61)
    --(75.65, 330.77)
    --(75.65, 327.50)
    --(75.65, 323.84)
    --(75.65, 319.81)
    --(75.65, 315.59)
    --(75.65, 311.24)
    --(75.65, 306.71)
    --(75.65, 302.24)
    --(75.65, 297.89)
    --(75.65, 293.67)
    --(75.65, 289.64)
    --(75.65, 285.98)
    --(75.65, 282.71)
    --(75.65, 279.87)
    --(75.65, 277.42)
    --(75.65, 275.40)
    --(75.65, 273.76)
    --(75.65, 272.57)
    --(75.65, 271.69)
    --(75.65, 271.18)
    --(75.65, 270.87)
    --(75.65, 270.74)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.75,0.75,0.75) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (80.13, 338.20)
    --(80.19, 338.20)
    --(80.44, 338.20)
    --(80.76, 338.27)
    --(81.32, 338.33)
    --(82.08, 338.39)
    --(83.09, 338.46)
    --(84.22, 338.52)
    --(85.61, 338.65)
    --(87.18, 338.77)
    --(88.82, 338.90)
    --(90.65, 339.02)
    --(92.54, 339.15)
    --(94.55, 339.28)
    --(96.63, 339.40)
    --(98.83, 339.46)
    --(101.10, 339.59)
    --(103.56, 339.72)
    --(106.14, 339.78)
    --(108.85, 339.84)
    --(111.81, 339.91)
    --(114.90, 339.97)
    --(118.17, 340.03)
    --(121.51, 340.03)
    --(124.85, 340.03)
    --(128.06, 339.97)
    --(131.15, 339.91)
    --(133.98, 339.84)
    --(136.69, 339.78)
    --(139.15, 339.72)
    --(141.48, 339.59)
    --(143.69, 339.46)
    --(145.76, 339.40)
    --(147.72, 339.28)
    --(149.61, 339.15)
    --(151.37, 339.02)
    --(153.07, 338.90)
    --(154.65, 338.77)
    --(156.03, 338.65)
    --(157.35, 338.52)
    --(158.43, 338.46)
    --(159.31, 338.39)
    --(160.00, 338.33)
    --(160.57, 338.27)
    --(160.88, 338.20)
    --(161.07, 338.20)
    --(161.13, 338.20)
    --(161.20, 338.27)
    --(161.39, 338.46)
    --(161.83, 338.90)
    --(162.52, 339.59)
    --(163.40, 340.47)
    --(164.22, 341.35)
    --(164.91, 342.05)
    --(165.35, 342.49)
    --(165.54, 342.68)
    --(165.61, 342.74)
    --(165.48, 342.74)
    --(165.23, 342.74)
    --(164.72, 342.74)
    --(163.97, 342.74)
    --(162.90, 342.74)
    --(161.45, 342.74)
    --(159.56, 342.74)
    --(157.29, 342.74)
    --(154.58, 342.74)
    --(151.50, 342.74)
    --(147.97, 342.74)
    --(144.06, 342.74)
    --(139.78, 342.74)
    --(135.24, 342.74)
    --(130.52, 342.74)
    --(125.61, 342.74)
    --(120.63, 342.74)
    --(115.65, 342.74)
    --(110.74, 342.74)
    --(106.02, 342.74)
    --(101.48, 342.74)
    --(97.20, 342.74)
    --(93.29, 342.74)
    --(89.76, 342.74)
    --(86.68, 342.74)
    --(83.97, 342.74)
    --(81.70, 342.74)
    --(79.81, 342.74)
    --(78.36, 342.74)
    --(77.29, 342.74)
    --(76.54, 342.74)
    --(76.03, 342.74)
    --(75.78, 342.74)
    --(75.65, 342.74)
    --(75.72, 342.68)
    --(75.91, 342.49)
    --(76.35, 342.05)
    --(77.04, 341.35)
    --(77.92, 340.47)
    --(78.74, 339.59)
    --(79.43, 338.90)
    --(79.87, 338.46)
    --(80.06, 338.27)
    --(80.13, 338.20)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.50,0.50,0.50) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (521.32, 449.95)
    --(521.39, 450.02)
    --(521.76, 450.39)
    --(522.52, 451.15)
    --(523.21, 451.91)
    --(523.59, 452.28)
    --(523.65, 452.35)
    --(523.65, 452.54)
    --(523.65, 452.91)
    --(523.65, 453.61)
    --(523.65, 454.68)
    --(523.65, 456.19)
    --(523.65, 458.08)
    --(523.65, 460.35)
    --(523.65, 462.93)
    --(523.65, 465.76)
    --(523.65, 468.79)
    --(523.65, 471.75)
    --(523.65, 474.58)
    --(523.65, 477.17)
    --(523.65, 479.43)
    --(523.65, 481.32)
    --(523.65, 482.83)
    --(523.65, 483.91)
    --(523.65, 484.60)
    --(523.65, 484.98)
    --(523.65, 485.17)
    --(523.59, 485.23)
    --(523.21, 485.61)
    --(522.46, 486.36)
    --(521.76, 487.06)
    --(521.39, 487.43)
    --(521.32, 487.50)
    --(521.32, 487.37)
    --(521.32, 486.99)
    --(521.32, 486.36)
    --(521.32, 485.35)
    --(521.32, 484.03)
    --(521.32, 482.27)
    --(521.32, 480.13)
    --(521.32, 477.61)
    --(521.32, 474.83)
    --(521.32, 471.81)
    --(521.32, 468.72)
    --(521.32, 465.64)
    --(521.32, 462.61)
    --(521.32, 459.84)
    --(521.32, 457.32)
    --(521.32, 455.18)
    --(521.32, 453.42)
    --(521.32, 452.09)
    --(521.32, 451.09)
    --(521.32, 450.46)
    --(521.32, 450.08)
    --(521.32, 449.95)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.75,0.75,0.75) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (523.65, 452.35)
    --(523.78, 452.35)
    --(523.97, 452.35)
    --(524.41, 452.28)
    --(525.04, 452.22)
    --(525.80, 452.16)
    --(526.80, 452.09)
    --(527.94, 452.03)
    --(529.20, 451.91)
    --(530.52, 451.84)
    --(531.97, 451.78)
    --(533.48, 451.65)
    --(535.06, 451.59)
    --(536.82, 451.53)
    --(538.71, 451.46)
    --(540.79, 451.46)
    --(542.99, 451.40)
    --(545.32, 451.40)
    --(547.65, 451.40)
    --(549.80, 451.46)
    --(551.81, 451.46)
    --(553.64, 451.53)
    --(555.34, 451.59)
    --(556.85, 451.65)
    --(558.30, 451.78)
    --(559.62, 451.84)
    --(560.88, 451.91)
    --(562.02, 452.03)
    --(563.02, 452.09)
    --(563.97, 452.16)
    --(564.66, 452.22)
    --(565.23, 452.28)
    --(565.61, 452.35)
    --(565.80, 452.35)
    --(565.92, 452.35)
    --(565.98, 452.28)
    --(566.36, 451.91)
    --(567.12, 451.15)
    --(567.87, 450.39)
    --(568.25, 450.02)
    --(568.31, 449.95)
    --(568.19, 449.95)
    --(567.81, 449.95)
    --(567.18, 449.95)
    --(566.17, 449.95)
    --(564.79, 449.95)
    --(563.02, 449.95)
    --(560.82, 449.95)
    --(558.17, 449.95)
    --(555.21, 449.95)
    --(551.87, 449.95)
    --(548.41, 449.95)
    --(544.82, 449.95)
    --(541.23, 449.95)
    --(537.76, 449.95)
    --(534.43, 449.95)
    --(531.46, 449.95)
    --(528.82, 449.95)
    --(526.61, 449.95)
    --(524.85, 449.95)
    --(523.46, 449.95)
    --(522.46, 449.95)
    --(521.83, 449.95)
    --(521.45, 449.95)
    --(521.32, 449.95)
    --(521.39, 450.02)
    --(521.76, 450.39)
    --(522.52, 451.15)
    --(523.21, 451.91)
    --(523.59, 452.28)
    --(523.65, 452.35)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((1.00,1.00,1.00) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (568.31, 449.95)
    --(568.25, 450.02)
    --(567.87, 450.39)
    --(567.12, 451.15)
    --(566.36, 451.91)
    --(565.98, 452.28)
    --(565.92, 452.35)
    --(565.92, 452.54)
    --(565.92, 452.91)
    --(565.92, 453.61)
    --(565.92, 454.68)
    --(565.92, 456.19)
    --(565.92, 458.08)
    --(565.92, 460.35)
    --(565.92, 462.93)
    --(565.92, 465.76)
    --(565.92, 468.79)
    --(565.92, 471.75)
    --(565.92, 474.58)
    --(565.92, 477.17)
    --(565.92, 479.43)
    --(565.92, 481.32)
    --(565.92, 482.83)
    --(565.92, 483.91)
    --(565.92, 484.60)
    --(565.92, 484.98)
    --(565.92, 485.17)
    --(565.98, 485.23)
    --(566.36, 485.61)
    --(567.12, 486.36)
    --(567.87, 487.06)
    --(568.25, 487.43)
    --(568.31, 487.50)
    --(568.31, 487.37)
    --(568.31, 486.99)
    --(568.31, 486.36)
    --(568.31, 485.35)
    --(568.31, 484.03)
    --(568.31, 482.27)
    --(568.31, 480.13)
    --(568.31, 477.61)
    --(568.31, 474.83)
    --(568.31, 471.81)
    --(568.31, 468.72)
    --(568.31, 465.64)
    --(568.31, 462.61)
    --(568.31, 459.84)
    --(568.31, 457.32)
    --(568.31, 455.18)
    --(568.31, 453.42)
    --(568.31, 452.09)
    --(568.31, 451.09)
    --(568.31, 450.46)
    --(568.31, 450.08)
    --(568.31, 449.95)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.75,0.75,0.75) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (523.65, 485.17)
    --(523.78, 485.17)
    --(523.97, 485.17)
    --(524.41, 485.23)
    --(525.04, 485.29)
    --(525.80, 485.35)
    --(526.80, 485.42)
    --(527.94, 485.54)
    --(529.20, 485.61)
    --(530.52, 485.67)
    --(531.97, 485.80)
    --(533.48, 485.86)
    --(535.06, 485.92)
    --(536.82, 486.05)
    --(538.71, 486.11)
    --(540.79, 486.11)
    --(542.99, 486.17)
    --(545.32, 486.17)
    --(547.65, 486.17)
    --(549.80, 486.11)
    --(551.81, 486.11)
    --(553.64, 486.05)
    --(555.34, 485.92)
    --(556.85, 485.86)
    --(558.30, 485.80)
    --(559.62, 485.67)
    --(560.88, 485.61)
    --(562.02, 485.54)
    --(563.02, 485.42)
    --(563.97, 485.35)
    --(564.66, 485.29)
    --(565.23, 485.23)
    --(565.61, 485.17)
    --(565.80, 485.17)
    --(565.92, 485.17)
    --(565.98, 485.23)
    --(566.36, 485.61)
    --(567.12, 486.36)
    --(567.87, 487.06)
    --(568.25, 487.43)
    --(568.31, 487.50)
    --(568.19, 487.50)
    --(567.81, 487.50)
    --(567.18, 487.50)
    --(566.17, 487.50)
    --(564.79, 487.50)
    --(563.02, 487.50)
    --(560.82, 487.50)
    --(558.17, 487.50)
    --(555.21, 487.50)
    --(551.87, 487.50)
    --(548.41, 487.50)
    --(544.82, 487.50)
    --(541.23, 487.50)
    --(537.76, 487.50)
    --(534.43, 487.50)
    --(531.46, 487.50)
    --(528.82, 487.50)
    --(526.61, 487.50)
    --(524.85, 487.50)
    --(523.46, 487.50)
    --(522.46, 487.50)
    --(521.83, 487.50)
    --(521.45, 487.50)
    --(521.32, 487.50)
    --(521.39, 487.43)
    --(521.76, 487.06)
    --(522.52, 486.30)
    --(523.21, 485.61)
    --(523.59, 485.23)
    --(523.65, 485.17)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.50,0.50,0.50) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (516.66, 414.80)
    --(558.93, 414.80)
    --(558.93, 413.92)
    --(516.66, 413.92)
    --(516.66, 414.80)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.50,0.50,0.50) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (520.44, 421.80)
    --(516.66, 414.80)
    --(523.65, 414.80)
    --(524.66, 417.13)
    --(523.02, 417.13)
    --(523.97, 419.46)
    --(522.46, 419.46)
    --(523.65, 421.80)
    --(520.44, 421.80)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.75,0.75,0.75) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (558.93, 414.80)
    --(551.87, 414.80)
    --(551.87, 417.13)
    --(553.39, 417.13)
    --(553.13, 419.46)
    --(554.20, 419.46)
    --(554.52, 419.46)
    --(554.08, 421.80)
    --(556.98, 421.80)
    --(558.93, 414.80)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.75,0.75,0.75) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (572.98, 414.80)
    --(570.83, 414.80)
    --(568.31, 420.41)
    --(559.31, 420.41)
    --(558.87, 421.80)
    --(569.26, 421.80)
    --(572.98, 414.80)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.75,0.75,0.75) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (561.26, 414.80)
    --(572.98, 414.80)
    --(572.98, 413.98)
    --(561.26, 413.98)
    --(561.26, 414.80)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.50,0.50,0.50) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=2;
  pickup pencircle scaled 0.00;
  path p;
  p = (561.26, 414.74)
    --(561.20, 413.92)
    --(558.87, 421.10)
    --(558.87, 421.80)
    --(561.26, 414.74)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.50,0.50,0.50) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=1;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (245.80, 691.53);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=1;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (247.12, 691.53);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=1;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (248.44, 691.53);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=1;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (249.76, 691.53);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=1;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (251.09, 691.53);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=1;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (246.49, 690.20);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=1;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (247.81, 690.20);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=1;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (249.13, 690.20);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=1;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (250.39, 690.20);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.00;
  path p;
  p = (330.02, 705.07)
    --(303.06, 714.08)
    --(213.04, 714.08)
    --(240.06, 705.07)
    --(330.02, 705.07)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.56,0.56,0.56) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (262.36, 699.91)
    --(264.63, 699.91);
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,1.00,0.00) + -0.25white);
  draw p withcolor (0.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (280.38, 699.91)
    --(282.65, 699.91);
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,1.00,0.00) + -0.25white);
  draw p withcolor (0.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (298.39, 699.91)
    --(300.60, 699.91);
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,1.00,0.00) + -0.25white);
  draw p withcolor (0.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (262.36, 693.17)
    --(264.63, 693.17);
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,1.00,0.00) + -0.25white);
  draw p withcolor (0.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (280.38, 693.17)
    --(282.65, 693.17);
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,1.00,0.00) + -0.25white);
  draw p withcolor (0.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (298.39, 693.17)
    --(300.60, 693.17);
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,1.00,0.00) + -0.25white);
  draw p withcolor (0.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (242.33, 696.13)
    --(244.54, 696.13);
  path f;
  f = p--cycle;
  fill f withcolor ((1.00,1.00,0.00) + -0.25white);
  draw p withcolor (1.00,1.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.00;
  path p;
  p = (240.06, 705.07)
    --(213.04, 714.08)
    --(213.04, 696.13)
    --(240.06, 687.12)
    --(240.06, 705.07)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.25,0.27,0.25) + 0.00white);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (85.92, 331.15)
    --(92.35, 331.15)
    --(103.37, 317.42)
    --(88.69, 298.96)
    --(85.92, 298.96)
    --(98.77, 315.53)
    --(85.92, 331.15)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((1.00,1.00,1.00) + 0.00white);
  draw p withcolor (1.00,1.00,1.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (118.93, 298.96)
    --(112.57, 298.96)
    --(101.48, 312.76)
    --(116.16, 331.15)
    --(118.93, 331.15)
    --(106.08, 314.58)
    --(118.93, 298.96)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((1.00,1.00,1.00) + 0.00white);
  draw p withcolor (1.00,1.00,1.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (258.33, 697.01)
    --(258.33, 698.33)
    --(260.60, 698.33)
    --(260.60, 697.01)
    --(260.16, 697.01)
    --(260.16, 696.57)
    --(259.21, 696.57)
    --(259.21, 697.01)
    --(258.33, 697.01)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,0.00,0.00) + -0.25white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (261.92, 697.01)
    --(261.92, 698.33)
    --(264.63, 698.33)
    --(264.63, 697.01)
    --(263.75, 697.01)
    --(263.75, 696.57)
    --(262.80, 696.57)
    --(262.80, 697.01)
    --(261.92, 697.01)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,0.00,0.00) + -0.25white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (269.61, 697.01)
    --(269.61, 698.33)
    --(271.81, 698.33)
    --(271.81, 697.01)
    --(270.93, 697.01)
    --(270.93, 696.57)
    --(270.05, 696.57)
    --(270.05, 697.01)
    --(269.61, 697.01)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,0.00,0.00) + -0.25white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (265.51, 697.01)
    --(265.51, 698.33)
    --(268.22, 698.33)
    --(268.22, 697.01)
    --(267.34, 697.01)
    --(267.34, 696.57)
    --(266.46, 696.57)
    --(266.46, 697.01)
    --(265.51, 697.01)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,0.00,0.00) + -0.25white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (276.35, 697.01)
    --(276.35, 698.33)
    --(278.55, 698.33)
    --(278.55, 697.01)
    --(278.11, 697.01)
    --(278.11, 696.57)
    --(277.23, 696.57)
    --(277.23, 697.01)
    --(276.35, 697.01)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,0.00,0.00) + -0.25white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (279.94, 697.01)
    --(279.94, 698.33)
    --(282.65, 698.33)
    --(282.65, 697.01)
    --(281.70, 697.01)
    --(281.70, 696.57)
    --(280.82, 696.57)
    --(280.82, 697.01)
    --(279.94, 697.01)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,0.00,0.00) + -0.25white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (287.56, 697.01)
    --(287.56, 698.33)
    --(289.83, 698.33)
    --(289.83, 697.01)
    --(288.94, 697.01)
    --(288.94, 696.57)
    --(288.00, 696.57)
    --(288.00, 697.01)
    --(287.56, 697.01)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,0.00,0.00) + -0.25white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (283.53, 697.01)
    --(283.53, 698.33)
    --(286.24, 698.33)
    --(286.24, 697.01)
    --(285.35, 697.01)
    --(285.35, 696.57)
    --(284.41, 696.57)
    --(284.41, 697.01)
    --(283.53, 697.01)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,0.00,0.00) + -0.25white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (294.30, 697.01)
    --(294.30, 698.33)
    --(296.57, 698.33)
    --(296.57, 697.01)
    --(296.13, 697.01)
    --(296.13, 696.57)
    --(295.24, 696.57)
    --(295.24, 697.01)
    --(294.30, 697.01)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,0.00,0.00) + -0.25white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (297.95, 697.01)
    --(297.95, 698.33)
    --(300.60, 698.33)
    --(300.60, 697.01)
    --(299.72, 697.01)
    --(299.72, 696.57)
    --(298.83, 696.57)
    --(298.83, 697.01)
    --(297.95, 697.01)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,0.00,0.00) + -0.25white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (305.57, 697.01)
    --(305.57, 698.33)
    --(307.84, 698.33)
    --(307.84, 697.01)
    --(306.90, 697.01)
    --(306.90, 696.57)
    --(306.02, 696.57)
    --(306.02, 697.01)
    --(305.57, 697.01)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,0.00,0.00) + -0.25white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (301.54, 697.01)
    --(301.54, 698.33)
    --(304.25, 698.33)
    --(304.25, 697.01)
    --(303.31, 697.01)
    --(303.31, 696.57)
    --(302.43, 696.57)
    --(302.43, 697.01)
    --(301.54, 697.01)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,0.00,0.00) + -0.25white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (258.33, 690.27)
    --(258.33, 691.59)
    --(260.60, 691.59)
    --(260.60, 690.27)
    --(260.16, 690.27)
    --(260.16, 689.83)
    --(259.21, 689.83)
    --(259.21, 690.27)
    --(258.33, 690.27)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,0.00,0.00) + -0.25white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (261.92, 690.27)
    --(261.92, 691.59)
    --(264.63, 691.59)
    --(264.63, 690.27)
    --(263.75, 690.27)
    --(263.75, 689.83)
    --(262.80, 689.83)
    --(262.80, 690.27)
    --(261.92, 690.27)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,0.00,0.00) + -0.25white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (269.61, 690.27)
    --(269.61, 691.59)
    --(271.81, 691.59)
    --(271.81, 690.27)
    --(270.93, 690.27)
    --(270.93, 689.83)
    --(270.05, 689.83)
    --(270.05, 690.27)
    --(269.61, 690.27)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,0.00,0.00) + -0.25white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (265.51, 690.27)
    --(265.51, 691.59)
    --(268.22, 691.59)
    --(268.22, 690.27)
    --(267.34, 690.27)
    --(267.34, 689.83)
    --(266.46, 689.83)
    --(266.46, 690.27)
    --(265.51, 690.27)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,0.00,0.00) + -0.25white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (276.35, 690.27)
    --(276.35, 691.59)
    --(278.55, 691.59)
    --(278.55, 690.27)
    --(278.11, 690.27)
    --(278.11, 689.83)
    --(277.23, 689.83)
    --(277.23, 690.27)
    --(276.35, 690.27)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,0.00,0.00) + -0.25white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (279.94, 690.27)
    --(279.94, 691.59)
    --(282.65, 691.59)
    --(282.65, 690.27)
    --(281.70, 690.27)
    --(281.70, 689.83)
    --(280.82, 689.83)
    --(280.82, 690.27)
    --(279.94, 690.27)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,0.00,0.00) + -0.25white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (287.56, 690.27)
    --(287.56, 691.59)
    --(289.83, 691.59)
    --(289.83, 690.27)
    --(288.94, 690.27)
    --(288.94, 689.83)
    --(288.00, 689.83)
    --(288.00, 690.27)
    --(287.56, 690.27)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,0.00,0.00) + -0.25white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (283.53, 690.27)
    --(283.53, 691.59)
    --(286.24, 691.59)
    --(286.24, 690.27)
    --(285.35, 690.27)
    --(285.35, 689.83)
    --(284.41, 689.83)
    --(284.41, 690.27)
    --(283.53, 690.27)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,0.00,0.00) + -0.25white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (294.30, 690.27)
    --(294.30, 691.59)
    --(296.57, 691.59)
    --(296.57, 690.27)
    --(296.13, 690.27)
    --(296.13, 689.83)
    --(295.24, 689.83)
    --(295.24, 690.27)
    --(294.30, 690.27)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,0.00,0.00) + -0.25white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (297.95, 690.27)
    --(297.95, 691.59)
    --(300.60, 691.59)
    --(300.60, 690.27)
    --(299.72, 690.27)
    --(299.72, 689.83)
    --(298.83, 689.83)
    --(298.83, 690.27)
    --(297.95, 690.27)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,0.00,0.00) + -0.25white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (305.57, 690.27)
    --(305.57, 691.59)
    --(307.84, 691.59)
    --(307.84, 690.27)
    --(306.90, 690.27)
    --(306.90, 689.83)
    --(306.02, 689.83)
    --(306.02, 690.27)
    --(305.57, 690.27)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,0.00,0.00) + -0.25white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (301.54, 690.27)
    --(301.54, 691.59)
    --(304.25, 691.59)
    --(304.25, 690.27)
    --(303.31, 690.27)
    --(303.31, 689.83)
    --(302.43, 689.83)
    --(302.43, 690.27)
    --(301.54, 690.27)--cycle;
  path f;
  f = p--cycle;
  fill f withcolor ((0.00,0.00,0.00) + -0.25white);
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (330.02, 687.12)
    --(330.02, 705.07)
    --(240.06, 705.07)
    --(240.06, 687.12)
    --(330.02, 687.12)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin text object
  picture q;
  q=thelabel.urt("Directory" infont "phvr" scaled (4pt/fontsize defaultfont),(0,0));
  picture p;
  p = q rotated 0.00;
  label.urt(p,((318.80,698.33))+xpart (lrcorner q - llcorner q)*(-cosd 0.00,-sind 0.00)/2+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  q=thelabel.urt("Service" infont "phvr" scaled (4pt/fontsize defaultfont),(0,0));
  picture p;
  p = q rotated 0.00;
  label.urt(p,((318.80,692.47))+xpart (lrcorner q - llcorner q)*(-cosd 0.00,-sind 0.00)/2+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
endfig;
end
