# thanks to http://patorjk.com/software/taag-v1/!

psychonetrics_print_logo <- function(color = c("default","red","green","yellow","blue","magenta","cyan","white","silver")){
  
  version <- read.dcf(file=system.file("DESCRIPTION", package="psychonetrics"),
                      fields="Version")
  
  version_string <- paste0("Version: ",version)
  
  # logo <- c("                        _                      _        _          ",
  #           "                      | |                    | |      (_)         ",
  #           "  _ __  ___ _   _  ___| |__   ___  _ __   ___| |_ _ __ _  ___ ___ ",
  #           " | '_ \\/ __| | | |/ __| '_ \\ / _ \\| '_ \\ / _ \\ __| '__| |/ __/ __|",
  #           " | |_) \\__ \\ |_| | (__| | | | (_) | | | |  __/ |_| |  | | (__\\__ \\",
  #           " | .__/|___/\\__, |\\___|_| |_|\\___/|_| |_|\\___|\\__|_|  |_|\\___|___/",
  #           " | |         __/ |                                                ",
  #           " |_|        |___/                                                 ","\n"
  # )
  color <- match.arg(color)
  logo <- c("                       _                      _        _          ",
            "                      | |                    | |      (_)         ",
            "  _ __  ___ _   _  ___| |__   ___  _ __   ___| |_ _ __ _  ___ ___ ",
            " |    \\/ __| | | |/ __|  _ \\ /   \\|  _ \\ / _ \\ __|  __| |/ __/ __|",
            " |  O--------------------------O  | | | |  __/ |_| |  | | (__\\__ \\",
            " | .__/|___/\\__, |\\___|_| |_|\\___/|_| |_|\\___|\\__|_|  |_|\\___|___/",
            " | |         __/ |                                                ",
            paste0(" |_|        |___/  ",paste0(rep(" ",66 - 19 - nchar(version_string)),collapse=""),version_string),
            "\n"
  )
if (color == "default"){
  cat(paste0(logo,collapse="\n"))
} else {
  cat(do.call(color,list(paste0(logo,collapse="\n"))))
}
  
}
