using Sockets

server = Sockets.listen(8080)

while true
    sock = Sockets.accept(server)
    @async begin
        data = readline(sock)
        print("Got request:\n", data, "\n") 
        cmd = split(data, " ")[2][2:end] 
        println(sock, "\nHTTP/1.1 200 OK\nContent-Type: text/html\n")
        println(sock, string("<html><body>", cmd, "=", eval(Meta.parse(cmd)), "</body></html>"))
        close(sock)
    end
end     