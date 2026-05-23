To implement Payne-Hanek reduction, we have to exploit a specific mathematical property of modulo arithmetic: **we don't need the full product of $x \times \frac{2}{\pi}$.**

Because $x$ has finite precision $P$ but a potentially massive exponent $E$, the higher bits of $\frac{2}{\pi}$ will be multiplied by the lowest bits of $x$ to create integers that are multiples of 4. Since we only care about the result modulo 4, we can completely ignore those high bits.

The core of arbitrary-precision Payne-Hanek relies on two steps:

1. **Precomputing $\frac{2}{\pi}$:** $2/\pi$ must be cached to the maximum supported exponent, stored as an integer bit-array.
2. **Slicing:** Instead of doing an $O(P \times E)$ multiplication, slice out a window of exactly $2P + W$ bits from the cached $2/\pi$ and do an $O(P \times P)$ multiplication.

### 1. The Global $2/\pi$ Cache

To make extracting the slice instant $O(1)$ instead of $O(E)$, store $2/\pi$ constant as a `big.Int` so we can slice its underlying `[]big.Word` array directly.

NOTE: we should actually pre-gen $2/\pi$ as data and store the []big.Word directly.

```go
var (
	// twoOverPiMaxBits is the absolute maximum exponent your library will support for Sin/Cos.
	twoOverPiMaxBits uint = 1_000_000 
	
	// twoOverPiBits holds the fractional bits of 2/pi.
	// It is conceptually equal to (2/pi) * 2^twoOverPiMaxBits
	twoOverPiBits *big.Int 
)

// Call this once during your library's init()
func initPayneHanekCache() {
	// Generate 2/pi to twoOverPiMaxBits of precision using your standard method
	f := twoOverPi(twoOverPiMaxBits) 
	
	// Shift the decimal point to the end to make it an integer
	f.SetMantExp(f, f.MantExp(nil)+int(twoOverPiMaxBits))
	
	twoOverPiBits = new(big.Int)
	f.Int(twoOverPiBits)
}

```

### 2. The $O(P)$ Slice Function

This function isolates the exact bit-window of $\frac{2}{\pi}$ that interacts with the mantissa of $x$ to produce the modulo 4 value and the fractional remainder.

```go
import "math/bits"

// sliceTwoOverPi returns the minimal window of 2/pi needed to reduce a number
// with exponent E, precision P, and guard bits G.
func sliceTwoOverPi(E, P, G uint) *big.Float {
	// The highest bit of 2/pi we must keep to preserve the modulo 4 result.
	// Any bit higher than this produces a multiple of 4 when multiplied by x.
	kStart := int(E) - int(P) - 2
	if kStart < 1 {
		kStart = 1 // highest bit of 2/pi is 2^-1
	}
	
	// The lowest bit needed to preserve the target precision
	kEnd := int(E) + int(P) + int(G)
	if kEnd > int(twoOverPiMaxBits) {
		panic("reducePi2: exponent exceeds max cached 2/pi bits")
	}

	numBits := uint(kEnd - kStart + 1)
	shift := twoOverPiMaxBits - uint(kEnd)

	// O(1) array slicing to find the relevant words
	words := twoOverPiBits.Bits()
	wordSize := uint(bits.UintSize)

	startWord := shift / wordSize
	startBit := shift % wordSize
	endWord := (shift + numBits + wordSize - 1) / wordSize
	if endWord > uint(len(words)) {
		endWord = uint(len(words))
	}

	var slice *big.Int
	if startWord >= uint(len(words)) {
		slice = big.NewInt(0)
	} else {
		// Copy only the ~O(P) required words 
		sliceWords := make([]big.Word, endWord-startWord)
		copy(sliceWords, words[startWord:endWord])

		slice = new(big.Int).SetBits(sliceWords)
		slice.Rsh(slice, startBit) // align exactly

		// Mask off any bits above the window
		mask := new(big.Int).Lsh(big.NewInt(1), numBits)
		mask.Sub(mask, big.NewInt(1))
		slice.And(slice, mask)
	}

	f := new(big.Float).SetPrec(numBits).SetInt(slice)

	// Restore the mathematical magnitude of the slice
	if f.Sign() != 0 {
		f.SetMantExp(f, f.MantExp(nil)-kEnd)
	}

	return f
}

```

### 3. Integrating it into `reducePi2`

Now we intercept the calculation. If $E$ is large enough that expanding `workPrec` would be costly, we switch to Payne-Hanek. Otherwise, we stick to the standard reduction because the slicing overhead isn't worth it for small numbers.

```go
func reducePi2(z, x *big.Float) int {
	if x.Sign() == 0 {
		return 0
	}

	prec := z.Prec()
	xAbs := newFloat(prec).Set(x)
	if xAbs.Signbit() {
		panic("reducePi2: x must be >= 0")
	}

	xExp := xAbs.MantExp(nil)
	workPrec := prec + _W

	var t0 *big.Float

	// Threshold heuristic: use Payne-Hanek if the exponent inflates 
	// precision by more than ~128 bits.
	if xExp > int(prec)+128 {
		// t0 = xAbs * Y_slice
		ySlice := sliceTwoOverPi(uint(xExp), prec, _W)
		
		// Notice workPrec stays extremely small (prec + _W)
		t0 = newFloat(workPrec).Mul(xAbs, ySlice)
	} else {
		if xExp > 0 {
			workPrec += uint(xExp)
		}
		// Standard full-precision fallback for small inputs
		t0 = newFloat(workPrec).Mul(xAbs, twoOverPi(workPrec))
	}

	// --- Proceed with the exact fast modulo logic from the previous step ---
	E := t0.MantExp(nil)
	var quadrant int

	if E <= 0 {
		t0.SetFloat64(0)
		quadrant = 0
	} else {
		t0.SetMode(big.ToZero)
		t0.SetPrec(uint(E))

		if E <= 2 {
			q64, _ := t0.Int64()
			quadrant = int(q64)
		} else {
			t1 := newFloat(uint(E - 2)).SetMode(big.ToZero).Set(t0)
			rem := newFloat(workPrec).Sub(t0, t1)
			q64, _ := rem.Int64()
			quadrant = int(q64)
		}
		t0.SetPrec(workPrec)
	}

	xAbs.Neg(xAbs)
	rTmp := fma(newFloat(workPrec), t0, halfPi(workPrec), xAbs, newFloat(workPrec))
	rTmp.Neg(rTmp)

	z.Set(rTmp)
	return quadrant
}

```
