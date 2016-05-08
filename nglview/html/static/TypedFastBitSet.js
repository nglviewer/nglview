/**
 * FastBitSet.js : a fast bit set implementation in JavaScript.
 * (c) the authors
 * Licensed under the Apache License, Version 2.0.
 *
 * Speed-optimized BitSet implementation for modern browsers and JavaScript engines.
 *
 * A BitSet is an ideal data structure to implement a Set when values being stored are
 * reasonably small integers. It can be orders of magnitude faster than a generic set implementation.
 * The FastBitSet implementation optimizes for speed, leveraging commonly available features
 * like typed arrays.
 *
 * Simple usage :
 *  var b = new TypedFastBitSet();// initially empty
 *         // will throw exception if typed arrays are not supported
 *  b.add(1);// add the value "1"
 *  b.has(1); // check that the value is present! (will return true)
 *  b.add(2);
 *  console.log(""+b);// should display {1,2}
 *  b.add(10);
 *  b.array(); // would return [1,2,10]
 *
 *  var c = new TypedFastBitSet([1,2,3,10]); // create bitset initialized with values 1,2,3,10
 *  c.difference(b); // from c, remove elements that are in b
 *  var su = c.union_size(b);// compute the size of the union (bitsets are unchanged)
 * c.union(b); // c will contain all elements that are in c and b
 * var s1 = c.intersection_size(b);// compute the size of the intersection (bitsets are unchanged)
 * c.intersection(b); // c will only contain elements that are in both c and b
 * c = b.clone(); // create a (deep) copy of b and assign it to c.
 * c.equals(b); // check whether c and b are equal
 *
 *   See README.md file for a more complete description.
 *
 * You can install the library under node with the command line
 *   npm install fastbitset
 */
'use strict';

// you can provide an iterable
// an exception is thrown if typed arrays are not supported
// - added size argument, ASR
// - added flip argument, ASR
function TypedFastBitSet(size, flip) {
  this.count = 0 | 0;
  this.words = new Uint32Array(8);
  this.resize(size);
  if (flip) {
    this.flip_all();
  }
}

// Add the value (Set the bit at index to true)
TypedFastBitSet.prototype.add = function(index) {
  if ((this.count << 5) <= index) {
    this.resize(index);
  }
  this.words[index >>> 5] |= 1 << index ;
};

// Add the value (Set the bit at index to true)
// - unsafe because size is not checked, added by ASR
TypedFastBitSet.prototype.add_unsafe = function(index) {
  this.words[index >>> 5] |= 1 << index ;
};

// If the value was not in the set, add it, otherwise remove it (flip bit at index)
TypedFastBitSet.prototype.flip = function(index) {
  if ((this.count << 5) <= index) {
    this.resize(index);
  }
  this.words[index >>> 5] ^= 1 << index ;
};

// If the value was not in the set, add it, otherwise remove it (flip bit at index)
// - unsafe because size is not checked, added by ASR
TypedFastBitSet.prototype.flip_unsafe = function(index) {
  this.words[index >>> 5] ^= 1 << index ;
};

// Flip all bits, added by ASR
TypedFastBitSet.prototype.flip_all = function() {
  // FIXME sets values beyond this.length
  var count = this.count;
  var k = 0 | 0;
  for (; k + 7 < count; k += 8) {
    this.words[k    ] = ~this.words[k    ];
    this.words[k + 1] = ~this.words[k + 1];
    this.words[k + 2] = ~this.words[k + 2];
    this.words[k + 3] = ~this.words[k + 3];
    this.words[k + 4] = ~this.words[k + 4];
    this.words[k + 5] = ~this.words[k + 5];
    this.words[k + 6] = ~this.words[k + 6];
    this.words[k + 7] = ~this.words[k + 7];
  }
  for (; k < count; ++k) {
    this.words[k] = ~this.words[k];
  }
  return this;
};

// Set all bits to value, added by ASR
TypedFastBitSet.prototype.set_all = function(value) {
  if (this.length <= 0) return this;
  this.set_range( 0, this.length-1, value );
  return this;
};

// Set all bits in range to value, added by ASR
TypedFastBitSet.prototype.set_range = function(from, to, value) {
  // set complete words when applicable
  value = value ? 0xFFFFFFFF : 0x00000000;
  var wordStart = (from >>> 5) + 1;
  var wordEnd = (to >>> 5) - 1;
  var k = wordStart | 0;
  for (; k + 7 < wordEnd; k += 8) {
    this.words[k    ] = value;
    this.words[k + 1] = value;
    this.words[k + 2] = value;
    this.words[k + 3] = value;
    this.words[k + 4] = value;
    this.words[k + 5] = value;
    this.words[k + 6] = value;
    this.words[k + 7] = value;
  }
  for (; k < wordEnd; ++k) {
    this.words[k] = value;
  }
  // set parts of the range not spanning complete words
  if (value){
    for (var i = from, n = (wordStart << 5) + 1; i < n; ++i) {
      this.words[i >>> 5] |= 1 << i ;
    }
    for (var i = (wordEnd << 5), n = to + 1; i < n; ++i) {
      this.words[i >>> 5] |= 1 << i ;
    }
  }else{
    for (var i = from, n = (wordStart << 5) + 1; i < n; ++i) {
      this.words[i >>> 5] &= ~(1 << i);
    }
    for (var i = (wordEnd << 5), n = to + 1; i < n; ++i) {
      this.words[i >>> 5] &= ~(1 << i);
    }
  }
  return this;
};

// Remove all values, reset memory usage
TypedFastBitSet.prototype.clear = function() {
  this.count = 0 | 0;
  this.length = 0 | 0;
  this.words = new Uint32Array(count);
  return this;
};

// Set the bit at index to false
TypedFastBitSet.prototype.remove = function(index) {
  if ((this.count << 5) <= index) {
    this.resize(index);
  }
  this.words[index >>> 5] &= ~(1 << index);
};

// Set the bit at index to false
// - unsafe because size is not checked, added by ASR
TypedFastBitSet.prototype.remove_unsafe = function(index) {
  this.words[index >>> 5] &= ~(1 << index);
};

// Return true if no bit is set
TypedFastBitSet.prototype.isEmpty = function(index) {
  var c = this.count;
  for (var  i = 0; i < c; i++) {
    if (this.words[i] !== 0) return false;
  }
  return true;
};

// Is the value contained in the set? Is the bit at index true or false? Returns a boolean
TypedFastBitSet.prototype.has = function(index) {
  return (this.words[index >>> 5] & (1 << index)) !== 0;
};

// Reduce the memory usage to a minimum
TypedFastBitSet.prototype.trim = function(index) {
  while (this.count > 0) {
    if (this.words[this.count - 1] === 0)
      this.count--;
  }
  this.words = this.words.slice(0,this.count);
};

// Resize the bitset so that we can write a value at index
TypedFastBitSet.prototype.resize = function(index) {
  this.length = index;
  if ((this.count << 5) > index) {
    return; //nothing to do
  }
  this.count = (index + 32) >>> 5;// just what is needed
  if ((this.words.length << 5) <= index) {
    var newwords = new Uint32Array(this.count << 1);
    newwords.set(this.words);// hopefully, this copy is fast
    this.words = newwords;
  }
};

// fast function to compute the Hamming weight of a 32-bit unsigned integer
TypedFastBitSet.prototype.hammingWeight = function(v) {
  v -= ((v >>> 1) & 0x55555555);// works with signed or unsigned shifts
  v = (v & 0x33333333) + ((v >>> 2) & 0x33333333);
  return ((v + (v >>> 4) & 0xF0F0F0F) * 0x1010101) >>> 24;
};


// How many values stored in the set? How many set bits?
TypedFastBitSet.prototype.size = function() {
  var answer = 0;
  var c = this.count;
  for (var i = 0; i < c; i++) {
    answer += this.hammingWeight(this.words[i] | 0);
  }
  return answer;
};

// How many bits are set in the given range of the set?
// - added by ASR
TypedFastBitSet.prototype.sizeRange = function(offset, count) {
  var size = 0;
  var end = offset + count;
  this.forEach( function( index ){
      if( index >= offset && index < end ) ++size;
  } );
  return size;
};

// Return an array with the set bit locations (values)
// - use Uint32Array instead of Array, ASR
TypedFastBitSet.prototype.array = function() {
  var answer = new Uint32Array(this.size());
  var pos = 0 | 0;
  var c = this.count | 0;
  for (var k = 0; k < c; ++k) {
    var w =  this.words[k];
    while (w != 0) {
      var t = w & -w;
      answer[pos++] = (k << 5) + this.hammingWeight((t - 1) | 0);
      w ^= t;
    }
  }
  return answer;
};


// Call fnc with the set bit locations (values)
// - fixed method description, ASR
TypedFastBitSet.prototype.forEach = function(fnc) {
  var c = this.count | 0;
  var i = 0 | 0;
  for (var k = 0; k < c; ++k) {
    var w = this.words[k];
    while (w != 0) {
      var t = w & -w;
      var index = (k << 5) + this.hammingWeight((t - 1) | 0);
      // FIXME workaround, it should not be required to check the length
      if(index < this.length) fnc(index,i);
      w ^= t;
      i += 1;
    }
  }
};

TypedFastBitSet.forEach = function(fnc, bitmap1, bitmap2) {
  var c = Math.min(bitmap1.count, bitmap2.count) | 0;
  for (var k = 0; k < c; ++k) {
    var w1 = bitmap1.words[k];
    var w2 = bitmap2.words[k];
    while (w1 != 0 && w2 != 0) {
      var t1 = w1 & -w1;
      var t2 = w2 & -w2;
      var kShift = k << 5;
      fnc(
        kShift + bitmap1.hammingWeight((t1 - 1) | 0),
        kShift + bitmap2.hammingWeight((t2 - 1) | 0)
      );
      w1 ^= t1;
      w2 ^= t2;
    }
  }
};

// Creates a copy of this bitmap
TypedFastBitSet.prototype.clone = function() {
  var clone = Object.create(TypedFastBitSet.prototype);
  clone.count = this.count;
  clone.length = this.length;
  clone.words = new Uint32Array(this.words);
  return clone;
};

// Check if this bitset intersects with another one,
// no bitmap is modified
TypedFastBitSet.prototype.intersects = function(otherbitmap) {
  var newcount = Math.min(this.count,otherbitmap.count);
  for (var k = 0 | 0; k < newcount; ++k) {
    if ((this.words[k] & otherbitmap.words[k]) !== 0) return true;
  }
  return false;
};

// Computes the intersection between this bitset and another one,
// the current bitmap is modified  (and returned by the function)
TypedFastBitSet.prototype.intersection = function(otherbitmap) {
  var newcount = Math.min(this.count,otherbitmap.count);
  var k = 0 | 0;
  for (; k + 7 < newcount; k += 8) {
    this.words[k    ] &= otherbitmap.words[k    ];
    this.words[k + 1] &= otherbitmap.words[k + 1];
    this.words[k + 2] &= otherbitmap.words[k + 2];
    this.words[k + 3] &= otherbitmap.words[k + 3];
    this.words[k + 4] &= otherbitmap.words[k + 4];
    this.words[k + 5] &= otherbitmap.words[k + 5];
    this.words[k + 6] &= otherbitmap.words[k + 6];
    this.words[k + 7] &= otherbitmap.words[k + 7];
  }
  for (; k < newcount; ++k) {
    this.words[k] &= otherbitmap.words[k];
  }
  var c = this.count;
  for (var k = newcount; k < c; ++k) {
    this.words[k] = 0;
  }
  this.count = newcount;
  return this;
};

// Computes the size of the intersection between this bitset and another one
TypedFastBitSet.prototype.intersection_size = function(otherbitmap) {
  var newcount = Math.min(this.count,otherbitmap.count);
  var answer = 0 | 0;
  for (var k = 0 | 0; k < newcount; ++k) {
    answer += this.hammingWeight(this.words[k] & otherbitmap.words[k]);
  }
  return answer;
};

// Computes the intersection between this bitset and another one,
// a new bitmap is generated
TypedFastBitSet.prototype.new_intersection = function(otherbitmap) {
  var answer = Object.create(TypedFastBitSet.prototype);
  answer.count = Math.min(this.count,otherbitmap.count);
  answer.words = new Uint32Array(answer.count);
  answer.length = Math.min(this.length,otherbitmap.length);
  var c = answer.count;
  for (var k = 0 | 0; k < c; ++k) {
    answer.words[k] = this.words[k] & otherbitmap.words[k];
  }
  return answer;
};

// Computes the intersection between this bitset and another one,
// the current bitmap is modified
TypedFastBitSet.prototype.equals = function(otherbitmap) {
  var mcount = Math.min(this.count , otherbitmap.count);
  for (var k = 0 | 0; k < mcount; ++k) {
    if (this.words[k] != otherbitmap.words[k]) return false;
  }
  if (this.count < otherbitmap.count) {
    var c = otherbitmap.count;
    for (var k = this.count; k < c; ++k) {
      if (otherbitmap.words[k] != 0) return false;
    }
  } else if (otherbitmap.count < this.count) {
    var c = this.count;
    for (var k = otherbitmap.count; k < c; ++k) {
      if (this.words[k] != 0) return false;
    }
  }
  return true;
};

// Computes the difference between this bitset and another one,
// the current bitset is modified (and returned by the function)
TypedFastBitSet.prototype.difference = function(otherbitmap) {
  var newcount = Math.min(this.count,otherbitmap.count);
  var k = 0 | 0;
  for (; k + 7 < newcount; k += 8) {
    this.words[k    ] &= ~otherbitmap.words[k    ];
    this.words[k + 1] &= ~otherbitmap.words[k + 1];
    this.words[k + 2] &= ~otherbitmap.words[k + 2];
    this.words[k + 3] &= ~otherbitmap.words[k + 3];
    this.words[k + 4] &= ~otherbitmap.words[k + 4];
    this.words[k + 5] &= ~otherbitmap.words[k + 5];
    this.words[k + 6] &= ~otherbitmap.words[k + 6];
    this.words[k + 7] &= ~otherbitmap.words[k + 7];
  }
  for (; k < newcount; ++k) {
    this.words[k] &= ~otherbitmap.words[k];
  }
  return this;
};

// Computes the size of the difference between this bitset and another one
TypedFastBitSet.prototype.difference_size = function(otherbitmap) {
  var newcount = Math.min(this.count,otherbitmap.count);
  var answer = 0 | 0;
  var k = 0 | 0;
  for (; k < newcount; ++k) {
    answer += this.hammingWeight(this.words[k] & (~otherbitmap.words[k]));
  }
  var c = this.count;
  for (; k < c; ++k) {
    answer += this.hammingWeight(this.words[k]);
  }
  return answer;
};

// Returns a string representation
TypedFastBitSet.prototype.toString = function() {
  return '{' + this.array().join(',') + '}';
};

// Computes the union between this bitset and another one,
// the current bitset is modified  (and returned by the function)
TypedFastBitSet.prototype.union = function(otherbitmap) {
  var mcount = Math.min(this.count,otherbitmap.count);
  var k = 0 | 0;
  for (; k + 7  < mcount; k += 8) {
    this.words[k    ] |= otherbitmap.words[k    ];
    this.words[k + 1] |= otherbitmap.words[k + 1];
    this.words[k + 2] |= otherbitmap.words[k + 2];
    this.words[k + 3] |= otherbitmap.words[k + 3];
    this.words[k + 4] |= otherbitmap.words[k + 4];
    this.words[k + 5] |= otherbitmap.words[k + 5];
    this.words[k + 6] |= otherbitmap.words[k + 6];
    this.words[k + 7] |= otherbitmap.words[k + 7];
  }
  for (; k < mcount; ++k) {
    this.words[k] |= otherbitmap.words[k];
  }
  if (this.count < otherbitmap.count) {
    this.resize((otherbitmap.count  << 5) - 1);
    var c = otherbitmap.count;
    for (var k = mcount; k < c; ++k) {
      this.words[k] = otherbitmap.words[k];
    }
    this.count = otherbitmap.count;
  }
  return this;
};

// Computes the union between this bitset and another one,
// a new bitmap is generated
TypedFastBitSet.prototype.new_union = function(otherbitmap) {
  var answer = Object.create(TypedFastBitSet.prototype);
  answer.count = Math.max(this.count,otherbitmap.count);
  answer.words = new Uint32Array(answer.count);
  var mcount = Math.min(this.count,otherbitmap.count)
  for (var k = 0; k < mcount; ++k) {
      answer.words[k] = this.words[k] | otherbitmap.words[k];
  }
  var c = this.count;
  for (var k = mcount; k < c; ++k) {
      answer.words[k] = this.words[k] ;
  }
  var c2 = otherbitmap.count;
  for (var k = mcount; k < c2; ++k) {
      answer.words[k] = otherbitmap.words[k] ;
  }
  return answer;
};

// Computes the difference between this bitset and another one,
// a new bitmap is generated
TypedFastBitSet.prototype.new_difference = function(otherbitmap) {
  return this.clone().difference(otherbitmap);// should be fast enough
};

// Computes the size union between this bitset and another one
TypedFastBitSet.prototype.union_size = function(otherbitmap) {
  var mcount = Math.min(this.count,otherbitmap.count);
  var answer = 0 | 0;
  for (var k = 0 | 0; k < mcount; ++k) {
    answer += this.hammingWeight(this.words[k] | otherbitmap.words[k]);
  }
  if (this.count < otherbitmap.count) {
    var c = otherbitmap.count;
    for (var k = this.count ; k < c; ++k) {
      answer += this.hammingWeight(otherbitmap.words[k] | 0);
    }
  } else {
    var c = this.count;
    for (var k = otherbitmap.count ; k < c; ++k) {
      answer += this.hammingWeight(this.words[k] | 0);
    }
  }
  return answer;
};

// Get transferable objects, added by ASR
TypedFastBitSet.prototype.getTransferable = function() {
  return [ this.words ];
};

// Serialize to JSON, added by ASR
TypedFastBitSet.prototype.toJSON = function() {
  return {
    count: this.count,
    length: this.length,
    words: this.words
  };
};

// De-serialize from JSON, added by ASR
TypedFastBitSet.prototype.fromJSON = function(input) {
  this.count = input.count;
  this.length = input.length;
  this.words = input.words;
  return this;
};
