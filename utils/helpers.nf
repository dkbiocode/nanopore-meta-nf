// utils/helpers.nf
c_bold = "\033[1m"
c_red = "\033[0;31m"
c_green = "\033[0;32m"
c_yellow = "\033[0;33m"
c_grey = "\033[90m"
c_reset = "\033[0m"

def checkDirs(Map dirMap, boolean createIfMissing = false) {
    dirMap.each { name, path ->
        if (file(path).exists()) {
            println("\t$name: $path;\t${c_green}exists${c_reset}")
        }
        else if (createIfMissing) {
            file(path).mkdirs()
            println("\t$name: $path;\t${c_yellow}created${c_reset}")
        }
        else {
            println("\t$name: $path;\t${c_red}MISSING${c_reset}")
        }
    }
}

def checkFiles(Map fileMap, boolean createDirIfMissing = false) {
    fileMap.each { name, path ->
        path_obj = file(path)
        if (path_obj.exists()) {
            println("\t$name:\t$path; ${c_green}exists${c_reset}")
        }
        else if (path_obj.parent.exists() ) {
            println("\t$name:\t${c_grey}$path${c_reset}; directory ${c_green}exists${c_reset} but file ${c_red}MISSING${c_reset}")
        }
        else if (createDirIfMissing) {
                path_obj.parent.mkdirs()
                println("\t$name:\t${c_grey}$path${c_reset}; directory ${c_yellow}created${c_reset} but file ${c_red}MISSING${c_reset}")
        }
        else {
            println("\t$name:\t$path; ${c_red}MISSING${c_reset}")
            println("\tcontaining dir: ${path_obj.parent};\t${c_red}MISSING${c_reset}")
        }
    }
}


