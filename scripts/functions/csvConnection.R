csvConnection <- function(x) {
    # source: https://stackoverflow.com/q/51448644
    # Open a connection to a raw vector
    write_rawcon <- rawConnection(object = raw(0), open = "w")

    # Write CSV file to the raw vector
    write.csv(x, file = write_rawcon, row.names = FALSE)

    # Access binary representation of CSV
    binary_csv <- rawConnectionValue(write_rawcon)

    # Close raw connection
    close(write_rawcon)

    # # Read binary back in
    # read_rawcon    <- rawConnection(object = binary_csv, open = "r")
    # new_binary_csv <- readBin(con = read_rawcon, what = "raw", n = length(binary_csv))

    # if (identical(binary_csv, new_binary_csv)) {
    #     return(new_binary_csv)
    # } else {
    #     stop("Something went wrong.")
    # }

    return(binary_csv)
}
