
#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(navbarPage("PALMS",
                   tabPanel("CalibrationCurve",
                            sidebarLayout(
                              sidebarPanel(
                                fileInput('file1', 'Choose CSV serial dilution table',
                                          accept=c('text/csv', 
                                                   'text/comma-separated-values,text/plain', 
                                                   '.csv')),
                                checkboxInput("gnps", "GnPS feature table", FALSE),
                                
                                conditionalPanel(
                                  condition = "input.gnps == true",
    				                      HTML("<font size=\"4\" color=\"red\">The order of the rows in the concentration table have to be in the same order as GnPS table columns.</font>")
                                ),
                                
                                conditionalPanel(
                                  condition = "input.gnps == true",
                                  fileInput('file11', 'Choose optional CSV concentration rownames',
                                            accept=c('text/csv', 
                                                     'text/comma-separated-values,text/plain', 
                                                     '.csv'))
                                ),
                                fileInput('file2', 'Choose optional CSV naming table',
                                          accept=c('text/csv', 
                                                   'text/comma-separated-values,text/plain', 
                                                   '.csv')),
                                checkboxInput("addpar", "Additional parameters.", FALSE),
    				                    conditionalPanel(
                                        condition = "input.addpar == true", 
            				                    numericInput("mz1", "m/z window (ppm)",
                                                     min = 1, max = 30, value = 10, step = 1),
                                        numericInput("rt1", "retention time window (seconds)",
                                                     min = 1, max = 60, value = 20, step = 1),
            				                    
            				                    
            				                    numericInput("slope", "Regression Slope threshold",
                                                     min = 0, max = 1, value = 0.1, step = 0.1),
            				                    
            				                    numericInput("r2", "Regression R2 threshold",
                                                     min = 0, max = 1, value = 0.9, step = 0.1),
            				                      
                                        checkboxInput("omass", "Allow only mass search", FALSE),
            				                      selectInput("ptype",
            				                                  "Plot Type:",
            				                                  c("regular", 
            				                                    "calplot")
            				                      )
    				                      ),
                                #actionButton("clustplot", "Plot cluster"),
                                br(), 
    				HTML("<font size=\"4\" color=\"black\">Use download for small feature tables only.</font>"), 
                                br(), 
                                downloadButton("data_file"),
                                
                                br(), 
                                br(), 
                                textInput("mail", "Email address:", ""),
                                actionButton("goMail", "Send calibration curve to email")
                                
                              ),
                              mainPanel(
                              )
                            )
                   ),
                   tabPanel("Quantification",
                            sidebarLayout(
                              sidebarPanel(
                                fileInput('file3', 'Choose .zip calibration curve',
                                          accept=c('.zip')),
                                fileInput('file4', 'Choose CSV feature table for quantification',
                                          accept=c('text/csv', 
                                                   'text/comma-separated-values,text/plain', 
                                                   '.csv')),
                                
                                numericInput("mz", "m/z window (ppm)",
                                             min = 1, max = 30, value = 10, step = 1),
                                numericInput("rt", "retention time window (seconds)",
                                             min = 1, max = 60, value = 20, step = 1),
                            
                                #actionButton("clustplot", "Plot cluster"),
                                br(), 
                                br(), 
    				HTML("<font size=\"4\" color=\"black\">Use download for small feature tables only.</font>"), 
                                br(), 
                                downloadButton("data_file2"),
				br(), 
                                br(), 
                                textInput("mail2", "Email address:", ""),
                                actionButton("goMail2", "Send quantitication table to email")
                              
                                
                              ),
                              mainPanel(
                              )
                            )
                   )
))
