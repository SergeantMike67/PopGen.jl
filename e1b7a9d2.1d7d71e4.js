(window.webpackJsonp=window.webpackJsonp||[]).push([[56],{156:function(e,t,a){"use strict";a.r(t),a.d(t,"frontMatter",(function(){return i})),a.d(t,"metadata",(function(){return s})),a.d(t,"rightToc",(function(){return c})),a.d(t,"default",(function(){return d}));var n=a(2),r=a(6),o=(a(0),a(165)),i={id:"datasets",title:"Provided datasets",sidebar_label:"Provided datasets"},s={id:"io/datasets",title:"Provided datasets",description:"PopGen.jl provides two datasets as examples, nancycats and gulfsharks. The datasets can be retrieved using the dataset function, or their names as commands without arguments (e.g. gulfsharks()).",source:"@site/docs/io/datasets.md",permalink:"/docs/io/datasets",editUrl:"https://github.com/pdimens/popgen.jl/edit/documentation/docs/io/datasets.md",sidebar_label:"Provided datasets",sidebar:"docs",previous:{title:"Variant Call Format",permalink:"/docs/io/vcf"},next:{title:"Start here",permalink:"/docs/tutorials/manipulate"}},c=[{value:"datasets",id:"datasets",children:[{value:"nancycats",id:"nancycats",children:[]},{value:"gulfsharks",id:"gulfsharks",children:[]}]}],l={rightToc:c};function d(e){var t=e.components,a=Object(r.a)(e,["components"]);return Object(o.b)("wrapper",Object(n.a)({},l,a,{components:t,mdxType:"MDXLayout"}),Object(o.b)("p",null,"PopGen.jl provides two datasets as examples, ",Object(o.b)("inlineCode",{parentName:"p"},"nancycats")," and ",Object(o.b)("inlineCode",{parentName:"p"},"gulfsharks"),". The datasets can be retrieved using the ",Object(o.b)("inlineCode",{parentName:"p"},"dataset")," function, or their names as commands without arguments (e.g. ",Object(o.b)("inlineCode",{parentName:"p"},"gulfsharks()"),"). "),Object(o.b)("div",{className:"admonition admonition-info alert alert--info"},Object(o.b)("div",Object(n.a)({parentName:"div"},{className:"admonition-heading"}),Object(o.b)("h5",{parentName:"div"},Object(o.b)("span",Object(n.a)({parentName:"h5"},{className:"admonition-icon"}),Object(o.b)("svg",Object(n.a)({parentName:"span"},{xmlns:"http://www.w3.org/2000/svg",width:"14",height:"16",viewBox:"0 0 14 16"}),Object(o.b)("path",Object(n.a)({parentName:"svg"},{fillRule:"evenodd",d:"M7 2.3c3.14 0 5.7 2.56 5.7 5.7s-2.56 5.7-5.7 5.7A5.71 5.71 0 0 1 1.3 8c0-3.14 2.56-5.7 5.7-5.7zM7 1C3.14 1 0 4.14 0 8s3.14 7 7 7 7-3.14 7-7-3.14-7-7-7zm1 3H6v5h2V4zm0 6H6v2h2v-2z"})))),"identitcal methods")),Object(o.b)("div",Object(n.a)({parentName:"div"},{className:"admonition-content"}),Object(o.b)("p",{parentName:"div"},"The methods are identical (one is a wrapper for the other), but the benefit of calling the datasets directly by name is that you get the luxury of tab auto-completion \ud83d\ude01"))),Object(o.b)("h2",{id:"datasets"},"datasets"),Object(o.b)("pre",null,Object(o.b)("code",Object(n.a)({parentName:"pre"},{className:"language-julia"}),"dataset(::String)\n")),Object(o.b)("p",null,"Returns a ",Object(o.b)("inlineCode",{parentName:"p"},"PopData")," object of the dataset you would like to retrieve by calling the dataset as a string by name."),Object(o.b)("p",null,Object(o.b)("strong",{parentName:"p"},"Example:")),Object(o.b)("pre",null,Object(o.b)("code",Object(n.a)({parentName:"pre"},{className:"language-julia"}),'sharks = dataset("gulfsharks")\ncats = dataset("nancycats")\n')),Object(o.b)("h3",{id:"nancycats"},"nancycats"),Object(o.b)("p",null,"We include the familiar nancycats microsatellite data, as featured in ",Object(o.b)("a",Object(n.a)({parentName:"p"},{href:"http://adegenet.r-forge.r-project.org"}),"adegenet"),", for easy importing into PopGen.jl as ",Object(o.b)("inlineCode",{parentName:"p"},"PopData"),". As an alternative to ",Object(o.b)("inlineCode",{parentName:"p"},"datasets"),", you can invoke the ",Object(o.b)("inlineCode",{parentName:"p"},"nancycats()"),"  command without any arguments."),Object(o.b)("pre",null,Object(o.b)("code",Object(n.a)({parentName:"pre"},{}),"julia> ncats = nancycats() ; summary(ncats)\nPopData Object\n  Marker type: Microsatellite\n  Ploidy: 2\n  Number of individuals: 237\n  Number of loci: 9\n  Populations: 17\n  Longitude: absent\n  Latitude: absent\n")),Object(o.b)("p",null,"The spatial coordinates provided for the dataset in ",Object(o.b)("inlineCode",{parentName:"p"},"adegenet")," are completely unfamiliar to us (and some geospatial folks we spoke to), so they have been omitted.  If you recognize what coordinate system has 485.111 appear in Nancy, France, please let us know!"),Object(o.b)("h3",{id:"gulfsharks"},"gulfsharks"),Object(o.b)("p",null,"We also include the SNP dataset used in Dimens ",Object(o.b)("em",{parentName:"p"},"et al."),' 2019 "',Object(o.b)("a",Object(n.a)({parentName:"p"},{href:"https://link.springer.com/article/10.1007/s00227-019-3533-1"}),"A genomic assessment of movement and gene flow around the South Florida vicariance zone in the migratory coastal blacknose shark, ",Object(o.b)("em",{parentName:"a"},"Carcharhinus acronotus")),'" since it was already on hand. Like ',Object(o.b)("inlineCode",{parentName:"p"},"nancycats"),", we provide a convenient function to load these data into PopGen.jl as ",Object(o.b)("inlineCode",{parentName:"p"},"PopData"),". As an alternative to ",Object(o.b)("inlineCode",{parentName:"p"},"datasets"),", you can invoke the ",Object(o.b)("inlineCode",{parentName:"p"},"gulfsharks()")," command without any arguments. "),Object(o.b)("pre",null,Object(o.b)("code",Object(n.a)({parentName:"pre"},{className:"language-julia"}),"julia> sharks = gulfsharks() ; summary(sharks)\nPopData Object\n  Marker type: SNP\n  Ploidy: 2\n  Number of individuals: 212\n  Number of loci: 2213\n  Populations: 7\n  Longitude: present with 0 missing\n  Latitude: present with 0 missing\n")))}d.isMDXComponent=!0},165:function(e,t,a){"use strict";a.d(t,"a",(function(){return p})),a.d(t,"b",(function(){return m}));var n=a(0),r=a.n(n);function o(e,t,a){return t in e?Object.defineProperty(e,t,{value:a,enumerable:!0,configurable:!0,writable:!0}):e[t]=a,e}function i(e,t){var a=Object.keys(e);if(Object.getOwnPropertySymbols){var n=Object.getOwnPropertySymbols(e);t&&(n=n.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),a.push.apply(a,n)}return a}function s(e){for(var t=1;t<arguments.length;t++){var a=null!=arguments[t]?arguments[t]:{};t%2?i(Object(a),!0).forEach((function(t){o(e,t,a[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(a)):i(Object(a)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(a,t))}))}return e}function c(e,t){if(null==e)return{};var a,n,r=function(e,t){if(null==e)return{};var a,n,r={},o=Object.keys(e);for(n=0;n<o.length;n++)a=o[n],t.indexOf(a)>=0||(r[a]=e[a]);return r}(e,t);if(Object.getOwnPropertySymbols){var o=Object.getOwnPropertySymbols(e);for(n=0;n<o.length;n++)a=o[n],t.indexOf(a)>=0||Object.prototype.propertyIsEnumerable.call(e,a)&&(r[a]=e[a])}return r}var l=r.a.createContext({}),d=function(e){var t=r.a.useContext(l),a=t;return e&&(a="function"==typeof e?e(t):s(s({},t),e)),a},p=function(e){var t=d(e.components);return r.a.createElement(l.Provider,{value:t},e.children)},b={inlineCode:"code",wrapper:function(e){var t=e.children;return r.a.createElement(r.a.Fragment,{},t)}},u=r.a.forwardRef((function(e,t){var a=e.components,n=e.mdxType,o=e.originalType,i=e.parentName,l=c(e,["components","mdxType","originalType","parentName"]),p=d(a),u=n,m=p["".concat(i,".").concat(u)]||p[u]||b[u]||o;return a?r.a.createElement(m,s(s({ref:t},l),{},{components:a})):r.a.createElement(m,s({ref:t},l))}));function m(e,t){var a=arguments,n=t&&t.mdxType;if("string"==typeof e||n){var o=a.length,i=new Array(o);i[0]=u;var s={};for(var c in t)hasOwnProperty.call(t,c)&&(s[c]=t[c]);s.originalType=e,s.mdxType="string"==typeof e?e:n,i[1]=s;for(var l=2;l<o;l++)i[l]=a[l];return r.a.createElement.apply(null,i)}return r.a.createElement.apply(null,a)}u.displayName="MDXCreateElement"}}]);