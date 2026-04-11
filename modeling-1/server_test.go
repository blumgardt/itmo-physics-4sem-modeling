package main

import (
	"bytes"
	"io"
	"net/http"
	"net/http/httptest"
	"strings"
	"testing"
)

func TestHandleComputeRejectsInvalidParams(t *testing.T) {
	req := httptest.NewRequest(http.MethodPost, "/api/compute", bytes.NewBufferString(`{"d":0.1,"l":0,"n":100,"i":1,"scale":1.5}`))
	rec := httptest.NewRecorder()

	handleCompute(rec, req)

	if rec.Code != http.StatusBadRequest {
		t.Fatalf("status = %d, want %d", rec.Code, http.StatusBadRequest)
	}

	body, err := io.ReadAll(rec.Body)
	if err != nil {
		t.Fatalf("read body: %v", err)
	}
	if !strings.Contains(string(body), "l must be > 0") {
		t.Fatalf("body = %q, want parameter validation error", string(body))
	}
}

func TestHandleComputeAcceptsMissingScale(t *testing.T) {
	req := httptest.NewRequest(http.MethodPost, "/api/compute", bytes.NewBufferString(`{"d":0.1,"l":0.3,"n":100,"i":1}`))
	rec := httptest.NewRecorder()

	handleCompute(rec, req)

	if rec.Code != http.StatusOK {
		t.Fatalf("status = %d, want %d", rec.Code, http.StatusOK)
	}
}
