def get_bits(number, lsb, offset):
    binary_representation = str(bin(number)[2:])[::-1]  # Convert the number to binary and remove '0b' prefix

    bits = (binary_representation[lsb:lsb+offset])[::-1]
    print("bianry rep: ", bits)
    # Extract the bits based on lsb and offset
    print("number: ",int(bits, 2) )
    return int(bits, 2)
  
def bits_to_represent(num):
    return  len(str(bin(num)[2:]))


def stupi_mul(k, P):
  window_size = 4 # window size
  n = bits_to_represent(k)
  
  print("n:", n)
  print("binary n: ", bin(k)[2:])
  
  precomp = [i*P for i in range(16)]
  
  print("precomp:", precomp, len(precomp))
  
  Q = 0
  m = int(n/window_size)
  
  for i in range(m):
    
    for j in range(window_size):
      Q = 2*Q

    precomp_idx = get_bits(k, int(m-i-1)*window_size, window_size)
    
    if (precomp_idx > 0):
      print("using precomp: ", precomp[precomp_idx])
      
      Q = Q + precomp[precomp_idx]
      
  return Q
  
print(stupi_mul(23029340, 23530))


