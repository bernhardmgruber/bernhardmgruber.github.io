---
layout: post
title:  "std::mdspan with struct of arrays (SoA) data layout"
date:   2023-11-22
hidden: true
---

I have been working on data layouts for a couple of years during my PhD.
Researching the relevant related work led me to `std::mdspan`,
or more accurately, [P0009](https://wg21.link/p0009), the C++ proposal for its adoption.
Given the customizable design of mdspan, which are mostly policies for layout and accessor,
I have been wondering for a while:

    Could mdspan provide a SoA layout?

Let's find out.

# Arrays of structs

Before digging into the matter, let's briefly talk about arrays of structs (AoS).
AoS is a data layout so ubiquitous we usually don't notice it.
Consider the following definitions to describe a pixel value:
```c++
struct RGBA {
   float r;
   float g
   float b;
   // 4 bytes padding
   double a;
};
```
When we now want to have an array of structs (see!), e.g. to describe an image, we can just say so:
```c++
struct Image {
    RGBA image[1024][1024];
};
```
Pixels and elements of the image are accesses using subscripts and data member access:
```
auto image = Image{};
auto& pixel = image[row][col];
pixel.r = 4.0;
```
And the data elements are laid out in memory as:
```
r, b, g, padding, a, r, g, b, padding, a, ...
```
That's it!
AoS is the "default" layout we get for arrays of structs (again!) in C++ and many other languages.

# The design of mdspan

Mdspan is a span.
As such, it does not own its memory (which is what mdarray will be for) but rather wraps an existing flat array of Ts,
which it represents as a multidimensional (MD) array, indexable by MD indices, forming an MD index space.
The shape of this MD index space is described by the space's extents.
Here is an example:

```c++
using Image = std::mdspan<std::extents<int, 1024, 1024>>;

std::vector<RGBA> storage(1024 * 1024);
auto image = Image{storage.data()};
auto& pixel = image[row, col]
pixel.r = 4.0;
```

Upon access, the passed MD index is mapped, i.e. translated, from the MD index space to a 1D index into the underlying flat array.
This aspect is customizable using the layout policy.
The standard provides some covering common use cases:

* `std::layout_right`: How C++ would lay out MD arrays, such as `int arr[3][4][5]`. For 2D arrays aka. row major.
* `std::layout_left`: How Fortran would lay out MD arrays, such as `integer, dimension (3, 4, 5) :: arr`. For 2D arrays aka. column major.
* `std::layout_stride`: Maps an MD index using `flat_index = dot_product(strides, md_index)`. This generalizes the upper two and can also skip elements. Useful to express slices of MD arrays. 
* and more! See e.g. [P2642](https://wg21.link/p2642)

After an MD index has been mapped to a 1D index, the underlying array can be subscripted.
By default, this is just a `underlying_array[flat_index]`, evaluating to a `T&`,
and this is in fact what `std::default_accessor<T>` does.
Here is a sketch:

```c++
template<typename T>
struct default_accessor {
  // ...
  using reference = T&;
  reference access(T* p, size_t i) const {
    return p[i];
  }
};
```

More accessors have been proposed, e.g. `std::atomic_accessor<T>` in [P2689](https://wg21.link/p2689),
which would return a `std::atomic_ref<T>` instead of T& as reference.

Because mdspan usually uses an underlying array of `T`s,
it will, independently of its layout policy, always use an array of structs (AoS) data layout.
There is some more flexibility given by the accessor's data handle,
which does not need to be a `T*`, but could be anything, e.g. a network handle to load data from.
However, this quickly complicates the setup of the underlying storage, so we will skip this for now. 

# Struct of arrays

Struct of arrays (SoA) is a data layout where the struct is split into its elements,
and those are laid out as separate arrays.
Reusing our image example, it could look like this:
```c++
struct Image {
   float r[1024][1024];
   float g[1024][1024];
   float b[1024][1024];
   double a[1024][1024];
};
```
Notice that there is no longer a struct to describe a single pixel (no `RGBA`).
And therefore, no reference to a pixel (no `RGBA&`).
Only references to elements can be created:
```c++
auto image = Image{};
auto& r = image.r[row][col];
// no construct to refer to pixel at [row][col] 
r = 4.0;
```
Also noticed, the indexing syntax swapped now.
We first select the pixel's element and then specify the array indices.

The data layout changed now, the elements are now laid out in memory as:
```c++
r, r, r, ..., g, g, g, ..., b, b, b, ..., a, a, a, ...
```
The contiguity of data elements of contiguity indices now allows for a bunch of optimizations,
notably vectorization, which is very important in number-crunching code.
Also notice, that no padding is required between data elements.
Also, all data is contained within one allocation.

There are many variations of SoA.
In practice, the array extents may not be known upfront, so each element may end up in its own allocation.
Our image may look something like this then:
```c++
struct Image {
   std::vector<float> r;
   std::vector<float> g;
   std::vector<float> b;
   std::vector<double> a;
};
```
There is a lot more about SoA, but we covered enough for now to follow this article.

# Struct of array reference

I mentioned that SoA layouts do not have an inherent reference type for their value type,
definitely not something as easy as `RGBA&` for an AoS layout.
A common trick is to use a struct of references for this purpose:
```c++
struct RGBARef {
    float& r;
    float& g;
    float& b;
    double& a;
};
```
Which we can create upon access using a little helper:

```c++
auto access(Image&, row, col) {
    return RGBARef{
        image.r[row][col],
        image.g[row][col],
        image.b[row][col],
        image.a[row][col]
    };
}

auto image = Image{};
auto pixel = access(image, row, col); // or: auto&&
pixel.r = 4.0;  // familiar syntax
```
We now got our familiar syntax back!
However, notice that we store the returned "reference" from the call to `access(...)` as value (`auto`),
and not as a `auto&`, which would be a compiler error (binding an l-value reference to a temporary).
A forwarding reference `auto&&` would also work, since it can hold onto the temporary **value** returned form `access(...)`.
You read that right, a SoA reference is a value in the type system, yet, modelling a reference.
Some people may call this a proxy reference, or a proxy object,
and this has all sorts of drawbacks (remember `std::vector<bool>`?)
But this is the price we pay for our familiar syntax.

We just scratched the surface here, since designing proper proxy reference types goes down a rabbit hole.
We may explore this in a different article :)

But why would we even want to have our familiar syntax?
(You can get out your C++ bingo cards now)
The answer is generic programming!
We would like to write a function that can work on both AoS and SoA layouts.
Therefore, the syntax has to look the same!
And mdspan already went ahead and unified the syntax for MD array access,
so we can write our matrix-manipulating functions without caring whether the matrix is stored row-major, column-major, triangular, sparse, etc.

# A SoA accessor

So, if mdspan gives us a uniform MD array syntax,
and our SoA proxy reference gives us a familiar member access syntax,
shouldn't we just be able to combine ... ? 
Yes, but before that, get ready for some UB :)
Since the underlying array of an mdspan is in AoS layout,
we have to completely write the memory addressing functionality ourselves and arbitrarily reinterpret bytes.

We can reuse the entire MD -> 1D layout mapping from mdspan,
so nothing to customize here.
Remember that the accessor had a nested `::reference` type, and an `access(...)` member function?
This is where we come in.
```c++
template<typename T, typename TRef>
struct SoAAccessor {
    using element_type = T;
    using reference = TRef;
```
An accessor also requires a nested `::data_handle_type` type,
which is usually a pointer to the underlying array (an AoS layout).
In our case, we completely disregard the structure of the underlying array and just abuse it as storage,
arbitrary reinterpreting the bytes.
```c++
    using data_handle_type = void*;
    static_assert(std::is_trivial_v<T>);
```
This is plain UB by the standard, even in cases where the static assertion passes.
But it works in practice for trivial types,
since they require no construction or destruction code.

And now we need a bunch of meta programming stuff,
for which I rely on [Boost.Mp11](https://www.boost.org/doc/libs/1_83_0/libs/mp11/doc/html/mp11.html).
Here are some helpers:
```c++
template <typename I>
constexpr auto roundUpToMultiple(I n, I mult) {
    return ((n + mult - 1) / mult) * mult;
}
```
This function rounds up an integral value to be a multiple of another one.

```c++
template<typename TypeList, std::size_t I>
inline constexpr std::size_t offsetOf = []() {
    if constexpr(I == 0)
        return 0;
    else {
        using T = boost::mp11::mp_at_c<TypeList, I - 1>;
        return roundUpToMultiple(offsetOf<TypeList, I - 1> + sizeof(T), alignof(T));
    }
}();
```
And this one is an equivalent of the [`offsetof`](https://en.cppreference.com/w/cpp/types/offsetof) macro,
but for type lists and using an index to select the member.
Example usage:
```c++
    using L = mp_list<float, float, float, double>; // RGBA struct as type list
    offsetOf<L, 0> ==  0;
    offsetOf<L, 1> ==  4;
    offsetOf<L, 2> ==  8;
    offsetOf<L, 3> == 16; // respects padding
```

To convert a struct into such a list, we can use [Boost.PFR](https://www.boost.org/doc/libs/master/doc/html/boost_pfr.html).
And we will place this line into our accessor as well:
```c++
    using Tuple = decltype(boost::pfr::structure_to_tuple(std::declval<T>()));
```

Furthermore, to compute our memory offsets,
we need to know the number of flat elements the mdspan stores.
This is given by `mdspan.mapping().required_span_size()`,
but we don't have access to the mapping from the accessor.
Now we may know that the mapping is just an adjacent data member of the accessor inside an mdspan instance ...
but maybe this is too much hacking :)
Furthermore, `required_span_size()` may involve computations to get the flat size from the extents,
which we don't want to have in our accessor's memory offset computation.
So let's just ask the user to provide it, and store it inside the accessor:
```c++
    int flatSize;
```

We are now ready to define the access function:
```c++
    reference access(data_handle_type p, std::size_t i) const noexcept {
        return [&]<std::size_t... Is>(std::index_sequence<Is...>) {
            return reference{
                *reinterpret_cast<std::tuple_element_t<Is, Tuple>*>(
                    reinterpret_cast<unsigned char*>(p)
                    + offsetOf<Tuple, Is> * flatSize
                    + sizeof(std::tuple_element_t<Is, Tuple>) * i
                )...
            };
        }(std::make_index_sequence<std::tuple_size_v<Tuple>>{});
    }
};
```
The immediately-invoked lambda expression just serves the purpose to create a pack of indices `Is` for the data members of our SoA reference type.
We can then construct our SoA reference `reference`, which is an aggregate btw.,
with a list of references to elements from the corresponding sub arrays.
Each such reference is computed as multiplying the data member offset by the flat number of elements,
which skips the preceding sub arrays,
adding the size of the element times the index,
and then adding all of that to the base pointer to our storage array.
All computations are done in bytes, so the storage array is reinterpreted as bytes,
and the final address is interpreted back as a pointer to the correct type.
This concludes our accessor implementation.

This implementation of a SoA data layout is by far not ideal.
It leaves a huge gap in the underlying storage where the padding was in the original `RGBA` struct:
```c++
r, r, r, ..., g, g, g, ..., b, b, b, ..., padding, padding, padding, ..., a, a, a, ...
```
It works though for the purpose of this article.
An improvement to avoid the padding is left as an exercise to the reader.
Hint: a space-efficient implementation is zero-overhead when `flatSize` is known at compile-time,
and probably requires a `roundUpToMultiple(...)` otherwise,
or the precomputation and storage of sub array offsets in the accesor instance.

# Evaluation

Let's construct an mdspan, but first with the `std::default_accessor`:
```c++
auto extents = std::dextents<int, 2>{1024, 1024};
auto mapping = std::layout_right::mapping{extents};
auto accessor = std::default_accessor<RGBA>{};
auto storage = std::vector<RGBA>(mapping.required_span_size());
auto image = std::mdspan{storage.data(), mapping, accessor};
```
We create extents, mapping, accessor, and a vector as underlying storage,
and then construct the mdspan on top.
Unfortunately, clang-17 failed to vectorize the code with static extents `std::extents<int, 1024, 1024>`,
so I am using dynamic extents here.

Now, we can call a small function to test our work.
It just scales the values of the red color channel:
```c++
void scaleRed(auto& image) {
    for (int row = 0; row < image.extent(1); row++)
        for (int col = 0; col < image.extent(0); col++) {
            auto&& pixel = image[row, col];
            pixel.r *= 1.5; // vector muls
        }
}

scaleRed(image);
```
Syntactically, it subscripts the mdspan, obtains a reference to a pixel, and accesses a data member of the pixel.

We compile using clang-17 and `-stdlib=libc++ -std=c++23 -O3 -mavx2`.
Inside the disassembly of `scaleRed(...)`, we can find the heart of the loop:
```asm
.LBB1_11:                               #   Parent Loop BB1_2 Depth=1
        vmovsd  xmm2, qword ptr [r12 - 32]      # xmm2 = mem[0],zero
        vmovhpd xmm2, xmm2, qword ptr [r12]     # xmm2 = xmm2[0],mem[0]
        vmovsd  xmm3, qword ptr [r12 - 96]      # xmm3 = mem[0],zero
        vmovhpd xmm3, xmm3, qword ptr [r12 - 64] # xmm3 = xmm3[0],mem[0]
        vinsertf128     ymm2, ymm3, xmm2, 1
        vmulpd  ymm2, ymm2, ymm1
        vmovlpd qword ptr [r12 - 96], xmm2
        vmovhpd qword ptr [r12 - 64], xmm2
        vextractf128    xmm2, ymm2, 1
        vmovlpd qword ptr [r12 - 32], xmm2
        vmovhpd qword ptr [r12], xmm2
        sub     r12, -128
        add     r13, -4
        jne     .LBB1_11
```
The compiler is already doing something very clever here.
It gathers the red elements from eight pixels, using a combination of instructions,
and puts them together into a single vector register (`ymm2`).
Then, it performs the multiplication to scale by `1.5` using `vmulpd`.
And then it scatters all the values back into memory.
At the end, it decrements `r13` by four, which is the loop counter,
restarting the loop if it has not hit zero yet.
The disassembly looks similar when compiling for AVX512 (using `-mavx512f`).
 
Now, let's change to our SoA accessor, reusing the exact same code otherwise:
```c++
auto accessor = AccessorSoA<RGBA, RGBARef>{mapping.required_span_size()};
```
Now, the disassembly of `scaleRed(...)` looks differently:
```asm
.LBB1_10:                               #   Parent Loop BB1_2 Depth=1
        vmulpd  ymm2, ymm1, ymmword ptr [r14 + 8*r15 - 96]
        vmulpd  ymm3, ymm1, ymmword ptr [r14 + 8*r15 - 64]
        vmulpd  ymm4, ymm1, ymmword ptr [r14 + 8*r15 - 32]
        vmulpd  ymm5, ymm1, ymmword ptr [r14 + 8*r15]
        vmovupd ymmword ptr [r14 + 8*r15 - 96], ymm2
        vmovupd ymmword ptr [r14 + 8*r15 - 64], ymm3
        vmovupd ymmword ptr [r14 + 8*r15 - 32], ymm4
        vmovupd ymmword ptr [r14 + 8*r15], ymm5
        add     r15, 16
        cmp     rbx, r15
        jne     .LBB1_10
```
And that looks like perfect vectorization.
Notice that no complicated gathering of red channel values is necessary.
The compiler, with each call to `vmulpd`, just loads 4 values and multiplies straight away.
The results are stored back to memory afterward,
with each `vmovupd` storing 4 values at the same time.
Notice also the four-fold duplication of SIMD instructions to allow for instruction-level parallelism.
This way, this loop pumps 16 doubles per iteration.
Using AVX512 (using `-mavx512f`), the compiler will emit the same code but using `zmm` registers,
processing 32 doubles per iteration.

Finally, I did a quick benchmark using Google benchmark on my AMD Ryzen with AVX2.
However, I took Kokkos' experimental mdspan implementation available via vcpkg together with libstdc++,
instead of the mdspan implementation in libc++ available since clang 17,
which I used on compiler explorer for the disassembly analysis.
The former was far easier to set up than convincing vcpkg to build Google benchmark with libc++.
It made a difference as we will see later.
```
-----------------------------------------------------
Benchmark           Time             CPU   Iterations
-----------------------------------------------------
AoS            803705 ns       803697 ns          781
SoA            110466 ns       110465 ns         6315
```
You can see that the SoA layout is about 8 times faster than AoS here.
I avoid a detailed analysis here for the sake of brevity,
but key contributors are likely:

* Better use of loaded memory.
  Each cacheline loaded from the AoS layout also contains data for the `g`, `b` and `a` channel, which we never need.
  The SoA does not load unnecessary data.
* No instructions wasted on reshuffling data (not the case on my machine, see below).
  The SoA loop runs pure compute, whereas the AoS has to do additional shuffling before it can crunch a single multiplication.
* Full use of vector registers for SoA instead of scalar code for AoS (my case).

As an aside, for the AoS version,
my local clang++ 17.0.4 using Kokkos' mdspan implementation even produced a scalar loop,
processing 2 doubles per iteration:
```asm
   vmulsd -0x20(%r11),%xmm0,%xmm1
   vmovsd %xmm1,-0x20(%r11)
   vmulsd (%r11),%xmm0,%xmm1
   vmovsd %xmm1,(%r11)
   add    $0x2,%r10
   add    $0x40,%r11
   cmp    %r10,%rsi
   jne    5ab10 <AoS(benchmark::State&)+0x110>
```

You can find the full code presented in this article on [compiler explorer](https://godbolt.org/z/vqvjo8jMG).
Additionally, I have used the following benchmark driver:
```c++
// insert code discussed so far, or from compiler explorer

static void AoS(benchmark::State& state) {
    auto extents = std::dextents<int, 2>{1024, 1024};
    auto mapping = std::layout_right::mapping{extents};
    auto accessor = std::default_accessor<RGBA>{};
    auto storage = std::vector<RGBA>(mapping.required_span_size());
    auto image = std::mdspan{storage.data(), mapping, accessor};

    for (auto _ : state) {
        scaleRed(image);
        benchmark::DoNotOptimize(image);
    }
}
BENCHMARK(AoS);

static void SoA(benchmark::State& state) {
    auto extents = std::dextents<int, 2>{1024, 1024};
    auto mapping = std::layout_right::mapping{extents};
    auto accessor = AccessorSoA<RGBA, RGBARef>{mapping.required_span_size()};
    auto storage = std::vector<RGBA>(mapping.required_span_size());
    auto image = std::mdspan{storage.data(), mapping, accessor};

    for (auto _ : state) {
        scaleRed(image);
        benchmark::DoNotOptimize(image);
    }
}
BENCHMARK(SoA);
```

# Another case for reflection

I want to end this article with a wish.
We desperately need proper reflection in C++.
Fortunately, [P2996](https://wg21.link/p2996r0) just dropped last month,
just in time for the 2023 WG21 fall meeting in Kona.
Even better: it was unanimously approved there by SG7, the study group for reflection,
and forwarded to EWG and LEWG.

If we had reflection:
* There would be no need for all the template metaprogramming arcanery to compute the offset of the i<sup>th</sup> data member of a struct.
  This would be a simple `std::meta::offset_of(std::meta::nonstatic_data_members_of(^T)[i])`.
* There is no need for the user to define a SoA reference, e.g. `RGBARef`, themselves,
  it could just be generated from the value type.
  This is analogous to the struct of arrays example in [P2996R0](https://wg21.link/p2996r0#struct-to-struct-of-arrays).

And the story does not end here, because my SoA reference sucks.
Imagine a pixel struct had operators defined for it,
or member functions,
which we suddenly would need to support on the distinct SoA reference type as well.
Reflection could just "copy" those operators to the new type.
Endless possibilities!

As a final word: Remember people saying that C++ gets more complicated with every release?
And how could it be otherwise, the committee just keeps adding things!
Well, reflection is one of those technologies that would just plainly eliminate a plethora of hacks we do these days,
making code indeed simpler, vastly simpler.
Let's hope it lands soon :)

# Acknowledgements

Many thanks to my friend Enrico for finding typos in this article.
Make sure to checkout his blog at [codekobold.io](https://codekobold.io/).
