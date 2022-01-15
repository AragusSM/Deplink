Below is a short list of the files and their purposes
This version of the app uses htmltemplate in conjunction with shiny.
The documentation is here. https://shiny.rstudio.com/articles/templates.html

Index.html:
The ui of the shiny server has been compiled into the index.html file, an html file
that also consists of R code in double brackets {{ "Rcode" }}
The index.html is contains many sections that represent the different parts of the webpage.
The system uses a hide show javascript system to hide the different div sections. 
The page is consisted of 5 div sections with the id indexdiv, aboutdiv, genediv, drugdiv, and contactdiv.
The page imports css and js files from online packages that are used to style and manipulate the page.
There are comment labels before each section to demarcate the pages and their content. 
Otherwise the code is fairly explanatory with html structure.

www: This folder contains the original logos and images. However, the assets folder contains
the majority of the javascript and css files. vendor contains core bootstrap and js files, and css
and js contains the custom designed js files for this version of the website and a package called
magnific popup (https://dimsemenov.com/plugins/magnific-popup/documentation.html#inline-type).
The images folder contains several stock images and some of the original images. Inside the team folder are
photos of staff members to put at the user's discretion. A sample image is provided. It is recommended to
resize additional photos to the sample dimensions.

server.ui:
Should be the same as before, should have all the shiny backend outputs and events.

ui.r:
simplified with the htmlTemplate() function. fill in the html parameters with shiny outputs and inputs.

data:
added a single image in tab 4 folder called message_plot2

function: 
unmodified

R files (rsconnect, .Rdata, etc):
should be unchanged as they are part of the project framework

