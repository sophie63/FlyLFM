# ----------------------------------------------------------------
#                 GRAYCODE DECODING UTILITIES
# ----------------------------------------------------------------

# Convert a binary code to a Gray code.
#
# This code from http://www.faqs.org/faqs/ai-faq/genetic/part6/section-1.html,
# pointed to from http://www.nist.gov/dads/HTML/graycode.html
def binaryToGraycode(binary_code, nbits):
    # Copy the high order bit
    gray_code = binary_code & (np.uint32(1) << nbits)
    
    # For the remaining bits, do the appropriate xor.
    shift = nbits-1
    while shift >= 0:
        Cb = (binary_code >> shift) & np.uint32(1)
        Cb1 = (binary_code >> (shift+1)) & np.uint32(1)
        gray_code |= ((Cb1^Cb) << shift)
        shift -= 1
        
    # The high order bits of each codeword are not used.
    return gray_code

def graycodeToBinary(gray_code, nbits):
    # Copy the high order bit
    binary_code = gray_code & (np.uint32(1) << nbits)
    
    # For the remaining bits, do the appropriate xor.
    shift = nbits-1
    while shift >= 0:
        
        Cg = (gray_code >> shift) & np.uint32(1)
        Cb1 = (binary_code >> (shift+1)) & np.uint32(1)
        binary_code |= ((Cb1^Cg) << shift)
        shift -= 1
        
    # The high order bits of each codeword are not used.
    return binary_code
