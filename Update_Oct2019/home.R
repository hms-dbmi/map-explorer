output$pageStub <- renderUI(tagList(
   fluidRow(
      column(10,
         HTML("<p>This is the home page of our excellent web site.</p>",
              "<p> You can upload your data here. </p>",
              "<p><p>There's a second page that displays summarized EHR information both on the overview level, 
              and on a detailed level with supporting evidence such as different lab measurement for the selected diseases for each patient.</p></p>",
              "<p> On that page you can ...</p>")
      )
   )
))
