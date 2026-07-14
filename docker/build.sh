docker buildx build \
  --platform linux/amd64,linux/arm64 \
  -t giacomoa/coding_challenge:1.0 \
  --push \
  . 
