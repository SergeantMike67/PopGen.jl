(window.webpackJsonp=window.webpackJsonp||[]).push([[20],{292:function(e,t,a){"use strict";a.r(t);var s=a(6),n=Object(s.a)({},(function(){var e=this,t=e.$createElement,a=e._self._c||t;return a("ContentSlotsDistributor",{attrs:{"slot-key":e.$parent.slotKey}},[a("h1",{attrs:{id:"other-data-types"}},[a("a",{staticClass:"header-anchor",attrs:{href:"#other-data-types"}},[e._v("#")]),e._v(" Other data types")]),e._v(" "),a("p",[e._v("While not strictly their own composite types, we also define aliases for genotypes and vectors of genotypes, as their explicit types can get a little unwieldy to use. The types shown below in the code blocks include their name and type (all types are of type "),a("code",[e._v("DataType")]),e._v(") on the first line, and what the alias actually defines on the second line.")]),e._v(" "),a("h3",{attrs:{id:"genotype"}},[a("a",{staticClass:"header-anchor",attrs:{href:"#genotype"}},[e._v("#")]),e._v(" Genotype")]),e._v(" "),a("div",{staticClass:"language-julia extra-class"},[a("pre",{pre:!0,attrs:{class:"language-julia"}},[a("code",[e._v("Genotype"),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v(":")]),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v(":")]),e._v("DataType\nNTuple"),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v("{")]),e._v("N"),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v(",")]),a("span",{pre:!0,attrs:{class:"token operator"}},[e._v("<:")]),e._v("Signed"),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v("}")]),e._v(" where N\n")])])]),a("p",[e._v("An "),a("code",[e._v("NTuple")]),e._v(" is itself an alias for a "),a("code",[e._v("Tuple{Vararg{}}")]),e._v(" , but you can think of it as Tuple of "),a("code",[e._v("N")]),e._v(" length composed of items of a particular type, in this case it's items that are subtypes of "),a("code",[e._v("Signed")]),e._v(" (the integer types). The length of the tuple ("),a("code",[e._v("N")]),e._v(") will vary based on the ploidy of the sample, and the element "),a("code",[e._v("Type")]),e._v(" will vary whether the markers are snps ("),a("code",[e._v("Int8")]),e._v(") or microsatellites ("),a("code",[e._v("Int16")]),e._v("), making this a pretty flexible (but immutable) structure.")]),e._v(" "),a("h3",{attrs:{id:"genotypearray"}},[a("a",{staticClass:"header-anchor",attrs:{href:"#genotypearray"}},[e._v("#")]),e._v(" GenotypeArray")]),e._v(" "),a("div",{staticClass:"language-julia extra-class"},[a("pre",{pre:!0,attrs:{class:"language-julia"}},[a("code",[e._v("GenotypeArray"),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v(":")]),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v(":")]),e._v("DataType\nAbstractVector"),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v("{")]),e._v("S"),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v("}")]),e._v(" where S"),a("span",{pre:!0,attrs:{class:"token operator"}},[e._v("<:")]),e._v("Union"),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v("{")]),e._v("Missing"),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v(",")]),e._v("Genotype"),a("span",{pre:!0,attrs:{class:"token punctuation"}},[e._v("}")]),e._v("\n")])])]),a("p",[e._v("As you can guess from the name, this Type wraps "),a("code",[e._v("Genotype")]),e._v(" into a Vector, while permitting "),a("code",[e._v("missing")]),e._v(" values (what's genetics without missing data!?). By using "),a("code",[e._v("AbstractVector")]),e._v(" (rather than "),a("code",[e._v("Vector")]),e._v("), we also have the flexibility of functions working on things like "),a("code",[e._v("SubArrays")]),e._v(" out of the box.")]),e._v(" "),a("div",{staticClass:"custom-block tip"},[a("p",{staticClass:"custom-block-title"},[e._v("why bother defining these aliases?")]),e._v(" "),a("p",[e._v("Getting the most out of Julia and demonstrating good practices means making sure functions work on the things they're supposed to, and give informative error messages when the input isn't suitable for the function (a rare case of "),a("em",[e._v("wanting")]),e._v(" MethodErrors). Without these aliases, functions would either have vague definitions like "),a("code",[e._v("f(x,y,z) where x <: AbstractArray")]),e._v(" and potentially cause errors, or overly complicated definitions like "),a("code",[e._v("f(x::AbstractVector{S},y,z) where {N, T<:Signed,S<:NTuple{N,T}}")]),e._v(" and not be very legible. Instead, functions are written as "),a("code",[e._v("f(x,y,z) where x<:GenotypeArray")]),e._v(", and that seems like a good compromise of getting the latter while looking like the former.")])])])}),[],!1,null,null,null);t.default=n.exports}}]);