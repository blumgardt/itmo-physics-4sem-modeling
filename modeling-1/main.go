package main

import (
	"flag"
	"fmt"
	"log"
	"net/http"
)

func main() {
	port := flag.Int("port", 8080, "HTTP server port")
	flag.Parse()

	http.HandleFunc("/", handleIndex)
	http.HandleFunc("/api/compute", handleCompute)

	addr := fmt.Sprintf(":%d", *port)
	fmt.Printf("Solenoid field simulator: http://localhost%s\n", addr)
	log.Fatal(http.ListenAndServe(addr, nil))
}
