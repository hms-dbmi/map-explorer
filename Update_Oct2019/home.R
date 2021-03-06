output$pageStub <- renderUI(tagList(
   fluidRow(
      column(10,
         h3('This is the home page of our excellent web app.'),
         h5('You can upload your data in this page.'),
         h5("Once you finish uploading all the required data, the original page will refresh and show you a new page that displays summarized EHR information both on the overview level, 
            and on a detailed level with supporting evidence such as different lab measurement for the selected diseases for each patient."),
         h5('On that page you can ...'),
         h5('You can always go back to the home page by clicking `Back to Home`'),
         
         h4("You can upload your own data files here!"),
         fileInput(inputId = "file1",
                   label = "Choose file for Phecode info:",
                   multiple = F,
                   buttonLabel = "Browse...",
                   placeholder = "No file selected"),
         fileInput(inputId = "file2"
                   , label = "Choose file for Patients and MAP cutoff:"
                   , multiple = F
                   , buttonLabel = "Browse..."
                   , placeholder = "No file selected"),
         fileInput(inputId = "file3"
                   , label = "Choose file for encounters of MS patient:"
                   , multiple = F
                   , buttonLabel = "Browse..."
                   , placeholder = "No file selected"),
         fileInput(inputId = "file4"
                   , label = "Choose file for Vitamin D values of MS patient:"
                   , multiple = F
                   , buttonLabel = "Browse..."
                   , placeholder = "No file selected"),
         fileInput(inputId = "file5"
                   , label = "Choose file for CUIs lists of MS patient:"
                   , multiple = F
                   , buttonLabel = "Browse..."
                   , placeholder = "No file selected")
      )
   )
))



