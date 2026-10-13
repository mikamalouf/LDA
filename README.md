# Luminex Data analysis (LDA)

Welcome to the Luminex data analysis app. This app is designed to help you analyse your Luminex data. Please upload your data using the file upload button on the left. You can upload multiple files at once. Once you have uploaded your data, you can use the tabs on the top of the app to visualise and process your data.

Please read before updating the code for this app.

The Shiny website can be found on https://jodyphelan.shinyapps.io/luminex/


# What does the app do?
1. Standard Curves
2. Variation <br />
  a. Coefficient of variation plots <br />
  b. Levey-Jennings plots <br />


**Definitions of each function** <br />
*ui* = user interface object controls the layout and appearance of your app <br />
*server* = function contains the instructions that your computer needs to build your app <br />
*page_sidebar* = function to create a page with a sidebar <br /> 
*fluidpage* = Functions for creating fluid page layouts. A fluid page layout consists of rows which in turn include columns <br />
*hr()* = Inserts a horizontal line (a divider in the UI) <br />
*h5("")* = adds level 5 heading <br />
*p()* = adds paragraph <br />

reactive({}) = R expression that uses widget input and returns a value. The reactive expression will update this value whenever the original widget changes.
- saves its result the first time you run it.
- The next time the reactive expression is called, it checks if the saved value has become out of date (i.e., whether the widgets it depends on have changed).
- If the value is out of date, the reactive object will recalculate it (and then save the new result).
- If the value is up-to-date, the reactive expression will return the saved value without doing any computation


# File format
1. Make sure that the date in the file is m/d/y <br />
2. The file upload must be a csv file <br />
3. Dilutions must be in order from highest to lowest concentration <br />
4. Dilutions must be in the format "1/100", not "100" or "1:100" or "0.01" <br />



# What to do before updating the code
1. Ensure that your R (different from RStudio) is up to date
2. Install required packages <br />

**Required packages** <br />
install.packages("coneproj") <br />
install.packages("cgam") <br />
install.packages("shinyFiles") <br />
install.packages("plotly") <br />
install.packages("plyr") <br />
install.packages("mixtools") <br />
install.packages("nls2") <br />
install.packages("lubridate")
