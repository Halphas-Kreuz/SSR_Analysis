#!/bin/bash

# Test for -gt (greater than)
if [ 5 -gt 3 ]; then
  echo "✅ -gt is working (5 is greater than 3)."
else
  echo "❌ -gt is NOT working."
fi

# Test for -eq (equal to)
if [ 5 -eq 5 ]; then
  echo "✅ -eq is working (5 is equal to 5)."
else
  echo "❌ -eq is NOT working."
fi