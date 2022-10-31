format_edits <- function(){
  shiny::tagList(
    shiny::tags$style(shiny::HTML(".skin-black .main-header .logo {background-color: white;}
                     .skin-black .main-header .logo:hover {background-color: white;}")),

    shiny::tags$style(shiny::HTML(".content-wrapper, .right-side{
      background-color: #002F6C;
      float: top;}


                                  ")),

    shiny::tags$head(shiny::tags$style(shiny::HTML(
      '.myClass {
        font-size: 20px;
        line-height: 50px;
        text-align: left;
        font-family: "Helvetica Neue",Helvetica,Arial,sans-serif;
        padding: 15px;
        overflow: hidden;
        color: #4D4D4D;
      }
    '
    )
    )
    ),

    shiny::tags$script(
      shiny::HTML('
        $(document).ready(function() {
        $("header").find("nav").append(\'<span class="myClass"> Survival Assessment Training and Practice</span>\');
      })
     ')
    )
  )}
