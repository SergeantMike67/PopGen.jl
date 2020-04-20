(window.webpackJsonp=window.webpackJsonp||[]).push([[15],{288:function(t,a,s){"use strict";s.r(a);var e=s(6),n=Object(e.a)({},(function(){var t=this,a=t.$createElement,s=t._self._c||a;return s("ContentSlotsDistributor",{attrs:{"slot-key":t.$parent.slotKey}},[s("h2",{attrs:{id:"delimited-format"}},[s("a",{staticClass:"header-anchor",attrs:{href:"#delimited-format"}},[t._v("#")]),t._v(" Delimited format")]),t._v(" "),s("h2",{attrs:{id:"import-a-delimited-file-as-popdata"}},[s("a",{staticClass:"header-anchor",attrs:{href:"#import-a-delimited-file-as-popdata"}},[t._v("#")]),t._v(" Import a delimited file as "),s("code",[t._v("PopData")])]),t._v(" "),s("div",{staticClass:"language-julia extra-class"},[s("pre",{pre:!0,attrs:{class:"language-julia"}},[s("code",[t._v("delimited"),s("span",{pre:!0,attrs:{class:"token punctuation"}},[t._v("(")]),t._v("infile"),s("span",{pre:!0,attrs:{class:"token punctuation"}},[t._v(":")]),s("span",{pre:!0,attrs:{class:"token punctuation"}},[t._v(":")]),t._v("String"),s("span",{pre:!0,attrs:{class:"token punctuation"}},[t._v(";")]),t._v(" delim"),s("span",{pre:!0,attrs:{class:"token punctuation"}},[t._v(":")]),s("span",{pre:!0,attrs:{class:"token punctuation"}},[t._v(":")]),t._v("Union"),s("span",{pre:!0,attrs:{class:"token punctuation"}},[t._v("{")]),t._v("Char"),s("span",{pre:!0,attrs:{class:"token punctuation"}},[t._v(",")]),t._v("String"),s("span",{pre:!0,attrs:{class:"token punctuation"}},[t._v(",")]),t._v("Regex"),s("span",{pre:!0,attrs:{class:"token punctuation"}},[t._v("}")]),t._v(" "),s("span",{pre:!0,attrs:{class:"token operator"}},[t._v("=")]),t._v(" "),s("span",{pre:!0,attrs:{class:"token string"}},[t._v('"auto"')]),s("span",{pre:!0,attrs:{class:"token punctuation"}},[t._v(",")]),t._v(" digits"),s("span",{pre:!0,attrs:{class:"token punctuation"}},[t._v(":")]),s("span",{pre:!0,attrs:{class:"token punctuation"}},[t._v(":")]),t._v("Int "),s("span",{pre:!0,attrs:{class:"token operator"}},[t._v("=")]),t._v(" "),s("span",{pre:!0,attrs:{class:"token number"}},[t._v("3")]),s("span",{pre:!0,attrs:{class:"token punctuation"}},[t._v(",")]),t._v(" diploid"),s("span",{pre:!0,attrs:{class:"token punctuation"}},[t._v(":")]),s("span",{pre:!0,attrs:{class:"token punctuation"}},[t._v(":")]),t._v("Bool "),s("span",{pre:!0,attrs:{class:"token operator"}},[t._v("=")]),t._v(" "),s("span",{pre:!0,attrs:{class:"token boolean"}},[t._v("true")]),s("span",{pre:!0,attrs:{class:"token punctuation"}},[t._v(",")]),t._v(" silent"),s("span",{pre:!0,attrs:{class:"token punctuation"}},[t._v(":")]),s("span",{pre:!0,attrs:{class:"token punctuation"}},[t._v(":")]),t._v("Bool "),s("span",{pre:!0,attrs:{class:"token operator"}},[t._v("=")]),t._v(" "),s("span",{pre:!0,attrs:{class:"token boolean"}},[t._v("false")]),s("span",{pre:!0,attrs:{class:"token punctuation"}},[t._v(")")]),t._v("\n\n"),s("span",{pre:!0,attrs:{class:"token comment"}},[t._v("# Example")]),t._v("\njulia"),s("span",{pre:!0,attrs:{class:"token operator"}},[t._v(">")]),t._v(" a "),s("span",{pre:!0,attrs:{class:"token operator"}},[t._v("=")]),t._v(" delimited"),s("span",{pre:!0,attrs:{class:"token punctuation"}},[t._v("(")]),s("span",{pre:!0,attrs:{class:"token string"}},[t._v('"/data/cali_poppy.csv"')]),s("span",{pre:!0,attrs:{class:"token punctuation"}},[t._v(",")]),t._v(" digits "),s("span",{pre:!0,attrs:{class:"token operator"}},[t._v("=")]),t._v(" "),s("span",{pre:!0,attrs:{class:"token number"}},[t._v("2")]),s("span",{pre:!0,attrs:{class:"token punctuation"}},[t._v(")")]),t._v("\n")])])]),s("div",{staticClass:"custom-block warning"},[s("p",{staticClass:"custom-block-title"},[t._v("Windows users")]),t._v(" "),s("p",[t._v("make sure to change your backslashes "),s("code",[t._v("\\")]),t._v(" to forward slashes "),s("code",[t._v("/")])])]),t._v(" "),s("h3",{attrs:{id:"arguments"}},[s("a",{staticClass:"header-anchor",attrs:{href:"#arguments"}},[t._v("#")]),t._v(" Arguments")]),t._v(" "),s("ul",[s("li",[s("code",[t._v("infile::String")]),t._v(" : path to the input file, in quotes")])]),t._v(" "),s("h3",{attrs:{id:"keyword-arguments"}},[s("a",{staticClass:"header-anchor",attrs:{href:"#keyword-arguments"}},[t._v("#")]),t._v(" Keyword Arguments")]),t._v(" "),s("ul",[s("li",[s("p",[s("code",[t._v("delim::String")]),t._v(" : delimiter characters. The default ("),s("code",[t._v('"auto"')]),t._v(") uses auto-parsing of "),s("code",[t._v("CSV.File")])])]),t._v(" "),s("li",[s("p",[s("code",[t._v("digits::Integer")]),t._v(" : the number of digits used to denote an allele (default: "),s("code",[t._v("3")]),t._v(")")])]),t._v(" "),s("li",[s("p",[s("code",[t._v("diploid::Bool")]),t._v("  : whether samples are diploid for parsing optimizations (default: "),s("code",[t._v("true")]),t._v(")")])]),t._v(" "),s("li",[s("p",[s("code",[t._v("silent::Bool")]),t._v(" : whether to print file information during import (default: "),s("code",[t._v("true")]),t._v(")")])])]),t._v(" "),s("h2",{attrs:{id:"formatting"}},[s("a",{staticClass:"header-anchor",attrs:{href:"#formatting"}},[t._v("#")]),t._v(" Formatting")]),t._v(" "),s("ul",[s("li",[t._v("First row is column names in this order:\n"),s("ol",[s("li",[t._v("name")]),t._v(" "),s("li",[t._v("population")]),t._v(" "),s("li",[t._v("longitude")]),t._v(" "),s("li",[t._v("latitude")]),t._v(" "),s("li",[t._v("locus_1_name")]),t._v(" "),s("li",[t._v("locus_2_name")]),t._v(" "),s("li",[t._v("etc...")])])])]),t._v(" "),s("h3",{attrs:{id:"missing-data"}},[s("a",{staticClass:"header-anchor",attrs:{href:"#missing-data"}},[t._v("#")]),t._v(" Missing data")]),t._v(" "),s("h4",{attrs:{id:"genotypes"}},[s("a",{staticClass:"header-anchor",attrs:{href:"#genotypes"}},[t._v("#")]),t._v(" Genotypes")]),t._v(" "),s("p",[t._v("Missing genotypes can be formatted as all-zeros (ex."),s("code",[t._v("000000")]),t._v(") or negative-nine "),s("code",[t._v("-9")])]),t._v(" "),s("h4",{attrs:{id:"location-data"}},[s("a",{staticClass:"header-anchor",attrs:{href:"#location-data"}},[t._v("#")]),t._v(" Location data")]),t._v(" "),s("p",[t._v("If location data is missing for a sample (which is ok!), make sure the value is written\nas "),s("code",[t._v("0")]),t._v(", otherwise there will be transcription errors!")]),t._v(" "),s("Tabs",{attrs:{card:"undefined",stretch:"undefined"}},[s("Tab",{attrs:{label:"formatting example"}},[s("div",{staticClass:"language- extra-class"},[s("pre",{pre:!0,attrs:{class:"language-text"}},[s("code",[t._v("name,population,long,lat,Locus1,Locus2,Locus3\nsierra_01,mountain,11.11,-22.22,001001,002002,001001\nsierra_02,mountain,11.12,-22.21,001001,001001,001002\nsnbarb_03,coast,0,0,001001,001001,001002\nsnbarb_02,coast,11.14,-22.24,001001,001001,001001\nsnbarb_03,coast,11.15,0,001002,001001,001001\n")])])])])],1),t._v(" "),s("div",{staticClass:"custom-block tip"},[s("p",{staticClass:"custom-block-title"},[t._v("TIP")]),t._v(" "),s("p",[t._v("You can also use the command "),s("code",[t._v("csv()")]),t._v(" synonymously with "),s("code",[t._v("delimited()")]),t._v(".")])]),t._v(" "),s("h2",{attrs:{id:"acknowledgements"}},[s("a",{staticClass:"header-anchor",attrs:{href:"#acknowledgements"}},[t._v("#")]),t._v(" Acknowledgements")]),t._v(" "),s("p",[t._v("Thanks to the efforts of the "),s("a",{attrs:{href:"https://github.com/JuliaData/CSV.jl",target:"_blank",rel:"noopener noreferrer"}},[t._v("CSV.jl"),s("OutboundLink")],1),t._v(" team, we are able leverage that package to do much of the heavy lifting within this parser.")])],1)}),[],!1,null,null,null);a.default=n.exports}}]);