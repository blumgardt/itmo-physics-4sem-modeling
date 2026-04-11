package main

import (
	_ "embed"
	"encoding/json"
	"fmt"
	"log"
	"math"
	"net/http"
)

//go:embed index.html
var indexHTML []byte

func validateParams(p Params) (Params, error) {
	if !isFinitePositive(p.D) {
		return p, fmt.Errorf("d must be > 0")
	}
	if !isFinitePositive(p.L) {
		return p, fmt.Errorf("l must be > 0")
	}
	if p.N < 1 {
		return p, fmt.Errorf("n must be >= 1")
	}
	if !isFinite(p.I) {
		return p, fmt.Errorf("i must be finite")
	}
	if p.Scale == 0 {
		p.Scale = 1.5
	}
	if !isFinite(p.Scale) || p.Scale < 1 {
		return p, fmt.Errorf("scale must be >= 1")
	}

	return p, nil
}

func isFinite(v float64) bool {
	return !math.IsNaN(v) && !math.IsInf(v, 0)
}

func isFinitePositive(v float64) bool {
	return isFinite(v) && v > 0
}

func handleIndex(w http.ResponseWriter, r *http.Request) {
	w.Header().Set("Content-Type", "text/html; charset=utf-8")
	w.Write(indexHTML)
}

func handleCompute(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodPost {
		http.Error(w, "POST only", http.StatusMethodNotAllowed)
		return
	}

	var p Params
	if err := json.NewDecoder(r.Body).Decode(&p); err != nil {
		http.Error(w, err.Error(), http.StatusBadRequest)
		return
	}
	p, err := validateParams(p)
	if err != nil {
		http.Error(w, err.Error(), http.StatusBadRequest)
		return
	}

	result := ComputeField(p)

	w.Header().Set("Content-Type", "application/json")
	if err := json.NewEncoder(w).Encode(result); err != nil {
		log.Printf("json encode error: %v", err)
	}
}
