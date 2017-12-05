runIoncopy <- 
function() {
  appDir <- system.file("shiny", "ioncopy", package="ioncopy")
  if (appDir == "") stop("Could not find ioncopy app directory!", call.=FALSE)
  shiny::runApp(appDir, display.mode="normal")
}
