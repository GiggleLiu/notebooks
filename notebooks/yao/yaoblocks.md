```@raw html
<style>
    table {
        display: table !important;
        margin: 2rem auto !important;
        border-top: 2pt solid rgba(0,0,0,0.2);
        border-bottom: 2pt solid rgba(0,0,0,0.2);
    }

    pre, div {
        margin-top: 1.4rem !important;
        margin-bottom: 1.4rem !important;
    }

    .code-output {
        padding: 0.7rem 0.5rem !important;
    }

    .admonition-body {
        padding: 0em 1.25em !important;
    }
</style>

<!-- PlutoStaticHTML.Begin -->
<!--
    # This information is used for caching.
    [PlutoStaticHTML.State]
    input_sha = "26d8c08c43eeac2da9630f7911f58d36381644c189868cad4953360bf501d2a1"
    julia_version = "1.8.0-rc1"
-->
<pre class='language-julia'><code class='language-julia'>TableOfContents()</code></pre>
<script>const getParentCell = el => el.closest("pluto-cell")

const getHeaders = () => {
	const depth = Math.max(1, Math.min(6, 3)) // should be in range 1:6
	const range = Array.from({length: depth}, (x, i) => i+1) // [1, ..., depth]
	
	const selector = range.map(i => `pluto-notebook pluto-cell h${i}`).join(",")
	return Array.from(document.querySelectorAll(selector))
}

const indent = true
const aside = true

const clickHandler = (event) => {
	const path = (event.path || event.composedPath())
	const toc = path.find(elem => elem?.classList?.contains?.("toc-toggle"))
	if (toc) {
		event.stopImmediatePropagation()
		toc.closest(".plutoui-toc").classList.toggle("hide")
	}
}

document.addEventListener("click", clickHandler)


const render = (el) => html`${el.map(h => {
	const parent_cell = getParentCell(h)

	const a = html`<a 
		class="${h.nodeName}" 
		href="#${parent_cell.id}"
	>${h.innerText}</a>`
	/* a.onmouseover=()=>{
		parent_cell.firstElementChild.classList.add(
			'highlight-pluto-cell-shoulder'
		)
	}
	a.onmouseout=() => {
		parent_cell.firstElementChild.classList.remove(
			'highlight-pluto-cell-shoulder'
		)
	} */
	a.onclick=(e) => {
		e.preventDefault();
		h.scrollIntoView({
			behavior: 'smooth', 
			block: 'start'
		})
	}

	return html`<div class="toc-row">${a}</div>`
})}`

const tocNode = html`<nav class="plutoui-toc">
	<header>
     <span class="toc-toggle open-toc">📖</span>
     <span class="toc-toggle closed-toc">📕</span>
	Table of Contents</header>
	<section></section>
</nav>`

tocNode.classList.toggle("aside", aside)
tocNode.classList.toggle("indent", indent)

const updateCallback = () => {
	tocNode.querySelector("section").replaceWith(
		html`<section>${render(getHeaders())}</section>`
	)
}
updateCallback()
setTimeout(updateCallback, 100)
setTimeout(updateCallback, 1000)
setTimeout(updateCallback, 5000)

const notebook = document.querySelector("pluto-notebook")


// We have a mutationobserver for each cell:
const observers = {
	current: [],
}

const createCellObservers = () => {
	observers.current.forEach((o) => o.disconnect())
	observers.current = Array.from(notebook.querySelectorAll("pluto-cell")).map(el => {
		const o = new MutationObserver(updateCallback)
		o.observe(el, {attributeFilter: ["class"]})
		return o
	})
}
createCellObservers()

// And one for the notebook's child list, which updates our cell observers:
const notebookObserver = new MutationObserver(() => {
	updateCallback()
	createCellObservers()
})
notebookObserver.observe(notebook, {childList: true})

// And finally, an observer for the document.body classList, to make sure that the toc also works when if is loaded during notebook initialization
const bodyClassObserver = new MutationObserver(updateCallback)
bodyClassObserver.observe(document.body, {attributeFilter: ["class"]})

// Hide/show the ToC when the screen gets small
let m = matchMedia("(max-width: 1000px)")
let match_listener = () => 
	tocNode.classList.toggle("hide", m.matches)
match_listener()
m.addListener(match_listener)

invalidation.then(() => {
	notebookObserver.disconnect()
	bodyClassObserver.disconnect()
	observers.current.forEach((o) => o.disconnect())
	document.removeEventListener("click", clickHandler)
	m.removeListener(match_listener)
})

return tocNode
</script><style>@media not print {

.plutoui-toc {
	--main-bg-color: unset;
	--pluto-output-color: hsl(0, 0%, 36%);
	--pluto-output-h-color: hsl(0, 0%, 21%);
}

@media (prefers-color-scheme: dark) {
	.plutoui-toc {
		--main-bg-color: hsl(0deg 0% 21%);
		--pluto-output-color: hsl(0, 0%, 90%);
		--pluto-output-h-color: hsl(0, 0%, 97%);
	}
}

.plutoui-toc.aside {
	color: var(--pluto-output-color);
	position:fixed;
	right: 1rem;
	top: 5rem;
	width: min(80vw, 300px);
	padding: 10px;
	border: 3px solid rgba(0, 0, 0, 0.15);
	border-radius: 10px;
	box-shadow: 0 0 11px 0px #00000010;
	/* That is, viewport minus top minus Live Docs */
	max-height: calc(100vh - 5rem - 56px);
	overflow: auto;
	z-index: 40;
	background-color: var(--main-bg-color);
	transition: transform 300ms cubic-bezier(0.18, 0.89, 0.45, 1.12);
}

.plutoui-toc.aside.hide {
	transform: translateX(calc(100% - 28px));
}
.plutoui-toc.aside.hide section {
	display: none;
}
.plutoui-toc.aside.hide header {
	margin-bottom: 0em;
	padding-bottom: 0em;
	border-bottom: none;
}
}  /* End of Media print query */
.plutoui-toc.aside.hide .open-toc,
.plutoui-toc.aside:not(.hide) .closed-toc,
.plutoui-toc:not(.aside) .closed-toc {
	display: none;
}

@media (prefers-reduced-motion) {
  .plutoui-toc.aside {
    transition-duration: 0s;
  }
}

.toc-toggle {
	cursor: pointer;
	padding: 1em;
	margin: -1em;
    margin-right: -0.7em;
}

.plutoui-toc header {
	display: block;
	font-size: 1.5em;
	margin-top: -0.1em;
	margin-bottom: 0.4em;
	padding-bottom: 0.4em;
	margin-left: 0;
	margin-right: 0;
	font-weight: bold;
	border-bottom: 2px solid rgba(0, 0, 0, 0.15);
}

.plutoui-toc section .toc-row {
	white-space: nowrap;
	overflow: hidden;
	text-overflow: ellipsis;
	padding-bottom: 2px;
}

.highlight-pluto-cell-shoulder {
	background: rgba(0, 0, 0, 0.05);
	background-clip: padding-box;
}

.plutoui-toc section a {
	text-decoration: none;
	font-weight: normal;
	color: var(--pluto-output-color);
}
.plutoui-toc section a:hover {
	color: var(--pluto-output-h-color);
}

.plutoui-toc.indent section a.H1 {
	font-weight: 700;
	line-height: 1em;
}

.plutoui-toc.indent section a.H1 {
	padding-left: 0px;
}
.plutoui-toc.indent section a.H2 {
	padding-left: 10px;
}
.plutoui-toc.indent section a.H3 {
	padding-left: 20px;
}
.plutoui-toc.indent section a.H4 {
	padding-left: 30px;
}
.plutoui-toc.indent section a.H5 {
	padding-left: 40px;
}
.plutoui-toc.indent section a.H6 {
	padding-left: 50px;
}
</style>

<pre class='language-julia'><code class='language-julia'>html"""
&lt;div align="center"&gt;
&lt;a class="Header-link " href="https://github.com/QuantumBFS/Yao.jl" data-hotkey="g d" aria-label="Homepage " data-ga-click="Header, go to dashboard, icon:logo"&gt;
  &lt;svg class="octicon octicon-mark-github v-align-middle" height="32" viewBox="0 0 16 16" version="1.1" width="32" aria-hidden="true"&gt;&lt;path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"&gt;&lt;/path&gt;&lt;/svg&gt;
&lt;/a&gt;
&lt;br&gt;
&lt;a href="https://raw.githubusercontent.com/GiggleLiu/notebooks/master/notebooks/yao/yaoblocks.jl" target="_blank"&gt; download this notebook &lt;/a&gt;&lt;/div&gt;
"""</code></pre>
<div align="center">
<a class="Header-link " href="https://github.com/QuantumBFS/Yao.jl" data-hotkey="g d" aria-label="Homepage " data-ga-click="Header, go to dashboard, icon:logo">
  <svg class="octicon octicon-mark-github v-align-middle" height="32" viewBox="0 0 16 16" version="1.1" width="32" aria-hidden="true"><path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"></path></svg>
</a>
<br>
<a href="https://raw.githubusercontent.com/GiggleLiu/notebooks/master/notebooks/yao/yaoblocks.jl" target="_blank"> download this notebook </a></div>



<div class="markdown"><h1>YaoBlocks - Hypercubic Linear algebra</h1>
</div>

<pre class='language-julia'><code class='language-julia'>using Yao, YaoPlots, PlutoUI, BenchmarkTools, KrylovKit</code></pre>



<div class="markdown"><h2>Definition and interfaces</h2>
<p>Unlike most quantum simulators, an Yao&#39;s block system is not only for a quantum gate, it implements the <strong>linear algebra on a hypercubic</strong>. All blocks subtypes <code>AbstractBlock</code>, and are classified into two categories</p>
<ol>
<li><p><code>PrimitiveBlock</code>: a leaf block that does not contain sub-blocks.</p>
</li>
<li><p><code>CompositeBlock</code>: a block composed of other blocks.</p>
</li>
</ol>
</div>

<pre class='language-julia'><code class='language-julia'>@with_terminal print_blocktree()</code></pre>
<pre id="plutouiterminal">
YaoAPI.AbstractBlock
    YaoAPI.CompositeBlock
        YaoAPI.AbstractContainer
            YaoAPI.TagBlock{BT} where BT<:YaoAPI.AbstractBlock
                YaoBlocks.AD.NoParams
                YaoBlocks.CachedBlock{ST, BT} where {ST, BT<:YaoAPI.AbstractBlock}
                YaoBlocks.Daggered
                YaoBlocks.Scale
                YaoPlots.LabelBlock
            YaoBlocks.ControlBlock
            YaoBlocks.PutBlock
            YaoBlocks.RepeatedBlock
            YaoBlocks.Subroutine
        YaoBlocks.AbstractAdd
            YaoBlocks.Add
        YaoBlocks.ChainBlock
        YaoBlocks.KronBlock
        YaoBlocks.UnitaryChannel
    YaoAPI.PrimitiveBlock
        Yao.EasyBuild.FSimGate
        YaoBlocks.ConstGate.ConstantGate
            Yao.EasyBuild.ISWAPGate
            Yao.EasyBuild.SqrtWGate
            Yao.EasyBuild.SqrtXGate
            Yao.EasyBuild.SqrtYGate
            YaoBlocks.ConstGate.CNOTGate
            YaoBlocks.ConstGate.CZGate
            YaoBlocks.ConstGate.HGate
            YaoBlocks.ConstGate.I2Gate
            YaoBlocks.ConstGate.IGate
            YaoBlocks.ConstGate.P0Gate
            YaoBlocks.ConstGate.P1Gate
            YaoBlocks.ConstGate.PdGate
            YaoBlocks.ConstGate.PuGate
            YaoBlocks.ConstGate.SGate
            YaoBlocks.ConstGate.SWAPGate
            YaoBlocks.ConstGate.SdagGate
            YaoBlocks.ConstGate.TGate
            YaoBlocks.ConstGate.TdagGate
            YaoBlocks.ConstGate.ToffoliGate
            YaoBlocks.ConstGate.XGate
            YaoBlocks.ConstGate.YGate
            YaoBlocks.ConstGate.ZGate
        YaoBlocks.GeneralMatrixBlock
        YaoBlocks.Measure
        YaoBlocks.PhaseGate
        YaoBlocks.RotationGate
        YaoBlocks.ShiftGate
        YaoBlocks.TimeEvolution
        YaoBlocks.TrivialGate
            YaoBlocks.IdentityGate
</pre>


<pre class='language-julia'><code class='language-julia'>qft3 = EasyBuild.qft_circuit(3)</code></pre>
<pre id='var-qft3' class='code-output documenter-example-output'>nqubits: 3
chain
├─ chain
│  ├─ put on (1)
│  │  └─ H
│  ├─ control(2)
│  │  └─ (1,) shift(1.5707963267948966)
│  └─ control(3)
│     └─ (1,) shift(0.7853981633974483)
├─ chain
│  ├─ put on (2)
│  │  └─ H
│  └─ control(3)
│     └─ (2,) shift(1.5707963267948966)
└─ chain
   └─ put on (3)
      └─ H
</pre>

<pre class='language-julia'><code class='language-julia'>vizcircuit(qft3)</code></pre>
<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:xlink="http://www.w3.org/1999/xlink"
     version="1.2"
     width="148.75mm" height="56.25mm" viewBox="0 0 148.75 56.25"
     stroke="none"
     fill="#000000"
     stroke-width="0.3"
     font-size="3.88"

     id="img-b91e1547">
<defs>
  <marker id="arrow" markerWidth="15" markerHeight="7" refX="5" refY="3.5" orient="auto" markerUnits="strokeWidth">
    <path d="M0,0 L15,3.5 L0,7 z" stroke="context-stroke" fill="context-stroke"/>
  </marker>
</defs>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-b91e1547-1">
  <g transform="translate(129.37,43.75)">
    <path d="M-5,-5 L5,-5 5,5 -5,5  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0" fill="#000000" id="img-b91e1547-2">
  <g transform="translate(108.12,43.75)">
    <circle cx="0" cy="0" r="0.83" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-b91e1547-3">
  <g transform="translate(108.12,25)">
    <path d="M-7.2,-5 L7.2,-5 7.2,5 -7.2,5  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-b91e1547-4">
  <g transform="translate(86.87,25)">
    <path d="M-5,-5 L5,-5 5,5 -5,5  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0" fill="#000000" id="img-b91e1547-5">
  <g transform="translate(65.62,43.75)">
    <circle cx="0" cy="0" r="0.83" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-b91e1547-6">
  <g transform="translate(65.62,6.25)">
    <path d="M-7.2,-5 L7.2,-5 7.2,5 -7.2,5  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0" fill="#000000" id="img-b91e1547-7">
  <g transform="translate(44.37,25)">
    <circle cx="0" cy="0" r="0.83" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-b91e1547-8">
  <g transform="translate(44.37,6.25)">
    <path d="M-7.2,-5 L7.2,-5 7.2,5 -7.2,5  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-b91e1547-9">
  <g transform="translate(23.12,6.25)">
    <path d="M-5,-5 L5,-5 5,5 -5,5  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-b91e1547-10">
  <g transform="translate(140,43.75)">
    <path fill="none" d="M0,0 L0,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-b91e1547-11">
  <g transform="translate(129.37,25)">
    <path fill="none" d="M10.62,0 L-10.62,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-b91e1547-12">
  <g transform="translate(108.12,6.25)">
    <path fill="none" d="M31.88,0 L-31.87,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-b91e1547-13">
  <g transform="translate(137.19,43.75)">
    <path fill="none" d="M-2.81,0 L2.81,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-b91e1547-14">
  <g transform="translate(121.56,43.75)">
    <path fill="none" d="M-2.81,0 L2.81,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-b91e1547-15">
  <g transform="translate(113.85,43.75)">
    <path fill="none" d="M-4.9,0 L4.9,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-b91e1547-16">
  <g transform="translate(91.77,43.75)">
    <path fill="none" d="M-15.52,0 L15.52,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-b91e1547-17">
  <g transform="translate(117.04,25)">
    <path fill="none" d="M-1.71,0 L1.71,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-b91e1547-18">
  <g transform="translate(99.21,25)">
    <path fill="none" d="M-1.71,0 L1.71,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-b91e1547-19">
  <g transform="translate(108.12,36.46)">
    <path fill="none" d="M0,6.46 L0,-6.46 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-b91e1547-20">
  <g transform="translate(94.69,25)">
    <path fill="none" d="M-2.81,0 L2.81,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-b91e1547-21">
  <g transform="translate(79.06,25)">
    <path fill="none" d="M-2.81,0 L2.81,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-b91e1547-22">
  <g transform="translate(65.62,25)">
    <path fill="none" d="M-10.62,0 L10.62,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-b91e1547-23">
  <g transform="translate(71.35,43.75)">
    <path fill="none" d="M-4.9,0 L4.9,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-b91e1547-24">
  <g transform="translate(38.65,43.75)">
    <path fill="none" d="M-26.15,0 L26.15,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-b91e1547-25">
  <g transform="translate(74.54,6.25)">
    <path fill="none" d="M-1.71,0 L1.71,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-b91e1547-26">
  <g transform="translate(56.71,6.25)">
    <path fill="none" d="M-1.71,0 L1.71,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-b91e1547-27">
  <g transform="translate(65.62,27.08)">
    <path fill="none" d="M0,15.83 L0,-15.83 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-b91e1547-28">
  <g transform="translate(50.1,25)">
    <path fill="none" d="M-4.9,0 L4.9,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-b91e1547-29">
  <g transform="translate(28.02,25)">
    <path fill="none" d="M-15.52,0 L15.52,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-b91e1547-30">
  <g transform="translate(53.29,6.25)">
    <path fill="none" d="M-1.71,0 L1.71,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-b91e1547-31">
  <g transform="translate(35.46,6.25)">
    <path fill="none" d="M-1.71,0 L1.71,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-b91e1547-32">
  <g transform="translate(44.37,17.71)">
    <path fill="none" d="M0,6.46 L0,-6.46 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-b91e1547-33">
  <g transform="translate(30.94,6.25)">
    <path fill="none" d="M-2.81,0 L2.81,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-b91e1547-34">
  <g transform="translate(15.31,6.25)">
    <path fill="none" d="M-2.81,0 L2.81,0 " class="primitive"/>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="5.64" id="img-b91e1547-35">
  <g transform="translate(129.37,43.75)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">H</text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="5.64" id="img-b91e1547-36">
  <g transform="translate(108.12,43.75)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em"></text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="3.53" id="img-b91e1547-37">
  <g transform="translate(108.12,25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">ϕ(π/2)</text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="5.64" id="img-b91e1547-38">
  <g transform="translate(86.87,25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">H</text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="5.64" id="img-b91e1547-39">
  <g transform="translate(65.62,43.75)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em"></text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="3.53" id="img-b91e1547-40">
  <g transform="translate(65.62,6.25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">ϕ(π/4)</text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="5.64" id="img-b91e1547-41">
  <g transform="translate(44.37,25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em"></text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="3.53" id="img-b91e1547-42">
  <g transform="translate(44.37,6.25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">ϕ(π/2)</text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="5.64" id="img-b91e1547-43">
  <g transform="translate(23.12,6.25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">H</text>
    </g>
  </g>
</g>
<script> <![CDATA[
(function(N){var k=/[\.\/]/,L=/\s*,\s*/,C=function(a,d){return a-d},a,v,y={n:{}},M=function(){for(var a=0,d=this.length;a<d;a++)if("undefined"!=typeof this[a])return this[a]},A=function(){for(var a=this.length;--a;)if("undefined"!=typeof this[a])return this[a]},w=function(k,d){k=String(k);var f=v,n=Array.prototype.slice.call(arguments,2),u=w.listeners(k),p=0,b,q=[],e={},l=[],r=a;l.firstDefined=M;l.lastDefined=A;a=k;for(var s=v=0,x=u.length;s<x;s++)"zIndex"in u[s]&&(q.push(u[s].zIndex),0>u[s].zIndex&&
(e[u[s].zIndex]=u[s]));for(q.sort(C);0>q[p];)if(b=e[q[p++] ],l.push(b.apply(d,n)),v)return v=f,l;for(s=0;s<x;s++)if(b=u[s],"zIndex"in b)if(b.zIndex==q[p]){l.push(b.apply(d,n));if(v)break;do if(p++,(b=e[q[p] ])&&l.push(b.apply(d,n)),v)break;while(b)}else e[b.zIndex]=b;else if(l.push(b.apply(d,n)),v)break;v=f;a=r;return l};w._events=y;w.listeners=function(a){a=a.split(k);var d=y,f,n,u,p,b,q,e,l=[d],r=[];u=0;for(p=a.length;u<p;u++){e=[];b=0;for(q=l.length;b<q;b++)for(d=l[b].n,f=[d[a[u] ],d["*"] ],n=2;n--;)if(d=
f[n])e.push(d),r=r.concat(d.f||[]);l=e}return r};w.on=function(a,d){a=String(a);if("function"!=typeof d)return function(){};for(var f=a.split(L),n=0,u=f.length;n<u;n++)(function(a){a=a.split(k);for(var b=y,f,e=0,l=a.length;e<l;e++)b=b.n,b=b.hasOwnProperty(a[e])&&b[a[e] ]||(b[a[e] ]={n:{}});b.f=b.f||[];e=0;for(l=b.f.length;e<l;e++)if(b.f[e]==d){f=!0;break}!f&&b.f.push(d)})(f[n]);return function(a){+a==+a&&(d.zIndex=+a)}};w.f=function(a){var d=[].slice.call(arguments,1);return function(){w.apply(null,
[a,null].concat(d).concat([].slice.call(arguments,0)))}};w.stop=function(){v=1};w.nt=function(k){return k?(new RegExp("(?:\\.|\\/|^)"+k+"(?:\\.|\\/|$)")).test(a):a};w.nts=function(){return a.split(k)};w.off=w.unbind=function(a,d){if(a){var f=a.split(L);if(1<f.length)for(var n=0,u=f.length;n<u;n++)w.off(f[n],d);else{for(var f=a.split(k),p,b,q,e,l=[y],n=0,u=f.length;n<u;n++)for(e=0;e<l.length;e+=q.length-2){q=[e,1];p=l[e].n;if("*"!=f[n])p[f[n] ]&&q.push(p[f[n] ]);else for(b in p)p.hasOwnProperty(b)&&
q.push(p[b]);l.splice.apply(l,q)}n=0;for(u=l.length;n<u;n++)for(p=l[n];p.n;){if(d){if(p.f){e=0;for(f=p.f.length;e<f;e++)if(p.f[e]==d){p.f.splice(e,1);break}!p.f.length&&delete p.f}for(b in p.n)if(p.n.hasOwnProperty(b)&&p.n[b].f){q=p.n[b].f;e=0;for(f=q.length;e<f;e++)if(q[e]==d){q.splice(e,1);break}!q.length&&delete p.n[b].f}}else for(b in delete p.f,p.n)p.n.hasOwnProperty(b)&&p.n[b].f&&delete p.n[b].f;p=p.n}}}else w._events=y={n:{}}};w.once=function(a,d){var f=function(){w.unbind(a,f);return d.apply(this,
arguments)};return w.on(a,f)};w.version="0.4.2";w.toString=function(){return"You are running Eve 0.4.2"};"undefined"!=typeof module&&module.exports?module.exports=w:"function"===typeof define&&define.amd?define("eve",[],function(){return w}):N.eve=w})(this);
(function(N,k){"function"===typeof define&&define.amd?define("Snap.svg",["eve"],function(L){return k(N,L)}):k(N,N.eve)})(this,function(N,k){var L=function(a){var k={},y=N.requestAnimationFrame||N.webkitRequestAnimationFrame||N.mozRequestAnimationFrame||N.oRequestAnimationFrame||N.msRequestAnimationFrame||function(a){setTimeout(a,16)},M=Array.isArray||function(a){return a instanceof Array||"[object Array]"==Object.prototype.toString.call(a)},A=0,w="M"+(+new Date).toString(36),z=function(a){if(null==
a)return this.s;var b=this.s-a;this.b+=this.dur*b;this.B+=this.dur*b;this.s=a},d=function(a){if(null==a)return this.spd;this.spd=a},f=function(a){if(null==a)return this.dur;this.s=this.s*a/this.dur;this.dur=a},n=function(){delete k[this.id];this.update();a("mina.stop."+this.id,this)},u=function(){this.pdif||(delete k[this.id],this.update(),this.pdif=this.get()-this.b)},p=function(){this.pdif&&(this.b=this.get()-this.pdif,delete this.pdif,k[this.id]=this)},b=function(){var a;if(M(this.start)){a=[];
for(var b=0,e=this.start.length;b<e;b++)a[b]=+this.start[b]+(this.end[b]-this.start[b])*this.easing(this.s)}else a=+this.start+(this.end-this.start)*this.easing(this.s);this.set(a)},q=function(){var l=0,b;for(b in k)if(k.hasOwnProperty(b)){var e=k[b],f=e.get();l++;e.s=(f-e.b)/(e.dur/e.spd);1<=e.s&&(delete k[b],e.s=1,l--,function(b){setTimeout(function(){a("mina.finish."+b.id,b)})}(e));e.update()}l&&y(q)},e=function(a,r,s,x,G,h,J){a={id:w+(A++).toString(36),start:a,end:r,b:s,s:0,dur:x-s,spd:1,get:G,
set:h,easing:J||e.linear,status:z,speed:d,duration:f,stop:n,pause:u,resume:p,update:b};k[a.id]=a;r=0;for(var K in k)if(k.hasOwnProperty(K)&&(r++,2==r))break;1==r&&y(q);return a};e.time=Date.now||function(){return+new Date};e.getById=function(a){return k[a]||null};e.linear=function(a){return a};e.easeout=function(a){return Math.pow(a,1.7)};e.easein=function(a){return Math.pow(a,0.48)};e.easeinout=function(a){if(1==a)return 1;if(0==a)return 0;var b=0.48-a/1.04,e=Math.sqrt(0.1734+b*b);a=e-b;a=Math.pow(Math.abs(a),
1/3)*(0>a?-1:1);b=-e-b;b=Math.pow(Math.abs(b),1/3)*(0>b?-1:1);a=a+b+0.5;return 3*(1-a)*a*a+a*a*a};e.backin=function(a){return 1==a?1:a*a*(2.70158*a-1.70158)};e.backout=function(a){if(0==a)return 0;a-=1;return a*a*(2.70158*a+1.70158)+1};e.elastic=function(a){return a==!!a?a:Math.pow(2,-10*a)*Math.sin(2*(a-0.075)*Math.PI/0.3)+1};e.bounce=function(a){a<1/2.75?a*=7.5625*a:a<2/2.75?(a-=1.5/2.75,a=7.5625*a*a+0.75):a<2.5/2.75?(a-=2.25/2.75,a=7.5625*a*a+0.9375):(a-=2.625/2.75,a=7.5625*a*a+0.984375);return a};
return N.mina=e}("undefined"==typeof k?function(){}:k),C=function(){function a(c,t){if(c){if(c.tagName)return x(c);if(y(c,"array")&&a.set)return a.set.apply(a,c);if(c instanceof e)return c;if(null==t)return c=G.doc.querySelector(c),x(c)}return new s(null==c?"100%":c,null==t?"100%":t)}function v(c,a){if(a){"#text"==c&&(c=G.doc.createTextNode(a.text||""));"string"==typeof c&&(c=v(c));if("string"==typeof a)return"xlink:"==a.substring(0,6)?c.getAttributeNS(m,a.substring(6)):"xml:"==a.substring(0,4)?c.getAttributeNS(la,
a.substring(4)):c.getAttribute(a);for(var da in a)if(a[h](da)){var b=J(a[da]);b?"xlink:"==da.substring(0,6)?c.setAttributeNS(m,da.substring(6),b):"xml:"==da.substring(0,4)?c.setAttributeNS(la,da.substring(4),b):c.setAttribute(da,b):c.removeAttribute(da)}}else c=G.doc.createElementNS(la,c);return c}function y(c,a){a=J.prototype.toLowerCase.call(a);return"finite"==a?isFinite(c):"array"==a&&(c instanceof Array||Array.isArray&&Array.isArray(c))?!0:"null"==a&&null===c||a==typeof c&&null!==c||"object"==
a&&c===Object(c)||$.call(c).slice(8,-1).toLowerCase()==a}function M(c){if("function"==typeof c||Object(c)!==c)return c;var a=new c.constructor,b;for(b in c)c[h](b)&&(a[b]=M(c[b]));return a}function A(c,a,b){function m(){var e=Array.prototype.slice.call(arguments,0),f=e.join("\u2400"),d=m.cache=m.cache||{},l=m.count=m.count||[];if(d[h](f)){a:for(var e=l,l=f,B=0,H=e.length;B<H;B++)if(e[B]===l){e.push(e.splice(B,1)[0]);break a}return b?b(d[f]):d[f]}1E3<=l.length&&delete d[l.shift()];l.push(f);d[f]=c.apply(a,
e);return b?b(d[f]):d[f]}return m}function w(c,a,b,m,e,f){return null==e?(c-=b,a-=m,c||a?(180*I.atan2(-a,-c)/C+540)%360:0):w(c,a,e,f)-w(b,m,e,f)}function z(c){return c%360*C/180}function d(c){var a=[];c=c.replace(/(?:^|\s)(\w+)\(([^)]+)\)/g,function(c,b,m){m=m.split(/\s*,\s*|\s+/);"rotate"==b&&1==m.length&&m.push(0,0);"scale"==b&&(2<m.length?m=m.slice(0,2):2==m.length&&m.push(0,0),1==m.length&&m.push(m[0],0,0));"skewX"==b?a.push(["m",1,0,I.tan(z(m[0])),1,0,0]):"skewY"==b?a.push(["m",1,I.tan(z(m[0])),
0,1,0,0]):a.push([b.charAt(0)].concat(m));return c});return a}function f(c,t){var b=O(c),m=new a.Matrix;if(b)for(var e=0,f=b.length;e<f;e++){var h=b[e],d=h.length,B=J(h[0]).toLowerCase(),H=h[0]!=B,l=H?m.invert():0,E;"t"==B&&2==d?m.translate(h[1],0):"t"==B&&3==d?H?(d=l.x(0,0),B=l.y(0,0),H=l.x(h[1],h[2]),l=l.y(h[1],h[2]),m.translate(H-d,l-B)):m.translate(h[1],h[2]):"r"==B?2==d?(E=E||t,m.rotate(h[1],E.x+E.width/2,E.y+E.height/2)):4==d&&(H?(H=l.x(h[2],h[3]),l=l.y(h[2],h[3]),m.rotate(h[1],H,l)):m.rotate(h[1],
h[2],h[3])):"s"==B?2==d||3==d?(E=E||t,m.scale(h[1],h[d-1],E.x+E.width/2,E.y+E.height/2)):4==d?H?(H=l.x(h[2],h[3]),l=l.y(h[2],h[3]),m.scale(h[1],h[1],H,l)):m.scale(h[1],h[1],h[2],h[3]):5==d&&(H?(H=l.x(h[3],h[4]),l=l.y(h[3],h[4]),m.scale(h[1],h[2],H,l)):m.scale(h[1],h[2],h[3],h[4])):"m"==B&&7==d&&m.add(h[1],h[2],h[3],h[4],h[5],h[6])}return m}function n(c,t){if(null==t){var m=!0;t="linearGradient"==c.type||"radialGradient"==c.type?c.node.getAttribute("gradientTransform"):"pattern"==c.type?c.node.getAttribute("patternTransform"):
c.node.getAttribute("transform");if(!t)return new a.Matrix;t=d(t)}else t=a._.rgTransform.test(t)?J(t).replace(/\.{3}|\u2026/g,c._.transform||aa):d(t),y(t,"array")&&(t=a.path?a.path.toString.call(t):J(t)),c._.transform=t;var b=f(t,c.getBBox(1));if(m)return b;c.matrix=b}function u(c){c=c.node.ownerSVGElement&&x(c.node.ownerSVGElement)||c.node.parentNode&&x(c.node.parentNode)||a.select("svg")||a(0,0);var t=c.select("defs"),t=null==t?!1:t.node;t||(t=r("defs",c.node).node);return t}function p(c){return c.node.ownerSVGElement&&
x(c.node.ownerSVGElement)||a.select("svg")}function b(c,a,m){function b(c){if(null==c)return aa;if(c==+c)return c;v(B,{width:c});try{return B.getBBox().width}catch(a){return 0}}function h(c){if(null==c)return aa;if(c==+c)return c;v(B,{height:c});try{return B.getBBox().height}catch(a){return 0}}function e(b,B){null==a?d[b]=B(c.attr(b)||0):b==a&&(d=B(null==m?c.attr(b)||0:m))}var f=p(c).node,d={},B=f.querySelector(".svg---mgr");B||(B=v("rect"),v(B,{x:-9E9,y:-9E9,width:10,height:10,"class":"svg---mgr",
fill:"none"}),f.appendChild(B));switch(c.type){case "rect":e("rx",b),e("ry",h);case "image":e("width",b),e("height",h);case "text":e("x",b);e("y",h);break;case "circle":e("cx",b);e("cy",h);e("r",b);break;case "ellipse":e("cx",b);e("cy",h);e("rx",b);e("ry",h);break;case "line":e("x1",b);e("x2",b);e("y1",h);e("y2",h);break;case "marker":e("refX",b);e("markerWidth",b);e("refY",h);e("markerHeight",h);break;case "radialGradient":e("fx",b);e("fy",h);break;case "tspan":e("dx",b);e("dy",h);break;default:e(a,
b)}f.removeChild(B);return d}function q(c){y(c,"array")||(c=Array.prototype.slice.call(arguments,0));for(var a=0,b=0,m=this.node;this[a];)delete this[a++];for(a=0;a<c.length;a++)"set"==c[a].type?c[a].forEach(function(c){m.appendChild(c.node)}):m.appendChild(c[a].node);for(var h=m.childNodes,a=0;a<h.length;a++)this[b++]=x(h[a]);return this}function e(c){if(c.snap in E)return E[c.snap];var a=this.id=V(),b;try{b=c.ownerSVGElement}catch(m){}this.node=c;b&&(this.paper=new s(b));this.type=c.tagName;this.anims=
{};this._={transform:[]};c.snap=a;E[a]=this;"g"==this.type&&(this.add=q);if(this.type in{g:1,mask:1,pattern:1})for(var e in s.prototype)s.prototype[h](e)&&(this[e]=s.prototype[e])}function l(c){this.node=c}function r(c,a){var b=v(c);a.appendChild(b);return x(b)}function s(c,a){var b,m,f,d=s.prototype;if(c&&"svg"==c.tagName){if(c.snap in E)return E[c.snap];var l=c.ownerDocument;b=new e(c);m=c.getElementsByTagName("desc")[0];f=c.getElementsByTagName("defs")[0];m||(m=v("desc"),m.appendChild(l.createTextNode("Created with Snap")),
b.node.appendChild(m));f||(f=v("defs"),b.node.appendChild(f));b.defs=f;for(var ca in d)d[h](ca)&&(b[ca]=d[ca]);b.paper=b.root=b}else b=r("svg",G.doc.body),v(b.node,{height:a,version:1.1,width:c,xmlns:la});return b}function x(c){return!c||c instanceof e||c instanceof l?c:c.tagName&&"svg"==c.tagName.toLowerCase()?new s(c):c.tagName&&"object"==c.tagName.toLowerCase()&&"image/svg+xml"==c.type?new s(c.contentDocument.getElementsByTagName("svg")[0]):new e(c)}a.version="0.3.0";a.toString=function(){return"Snap v"+
this.version};a._={};var G={win:N,doc:N.document};a._.glob=G;var h="hasOwnProperty",J=String,K=parseFloat,U=parseInt,I=Math,P=I.max,Q=I.min,Y=I.abs,C=I.PI,aa="",$=Object.prototype.toString,F=/^\s*((#[a-f\d]{6})|(#[a-f\d]{3})|rgba?\(\s*([\d\.]+%?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+%?(?:\s*,\s*[\d\.]+%?)?)\s*\)|hsba?\(\s*([\d\.]+(?:deg|\xb0|%)?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+(?:%?\s*,\s*[\d\.]+)?%?)\s*\)|hsla?\(\s*([\d\.]+(?:deg|\xb0|%)?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+(?:%?\s*,\s*[\d\.]+)?%?)\s*\))\s*$/i;a._.separator=
RegExp("[,\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]+");var S=RegExp("[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*"),X={hs:1,rg:1},W=RegExp("([a-z])[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029,]*((-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*)+)",
"ig"),ma=RegExp("([rstm])[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029,]*((-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*)+)","ig"),Z=RegExp("(-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?)[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*",
"ig"),na=0,ba="S"+(+new Date).toString(36),V=function(){return ba+(na++).toString(36)},m="http://www.w3.org/1999/xlink",la="http://www.w3.org/2000/svg",E={},ca=a.url=function(c){return"url('#"+c+"')"};a._.$=v;a._.id=V;a.format=function(){var c=/\{([^\}]+)\}/g,a=/(?:(?:^|\.)(.+?)(?=\[|\.|$|\()|\[('|")(.+?)\2\])(\(\))?/g,b=function(c,b,m){var h=m;b.replace(a,function(c,a,b,m,t){a=a||m;h&&(a in h&&(h=h[a]),"function"==typeof h&&t&&(h=h()))});return h=(null==h||h==m?c:h)+""};return function(a,m){return J(a).replace(c,
function(c,a){return b(c,a,m)})}}();a._.clone=M;a._.cacher=A;a.rad=z;a.deg=function(c){return 180*c/C%360};a.angle=w;a.is=y;a.snapTo=function(c,a,b){b=y(b,"finite")?b:10;if(y(c,"array"))for(var m=c.length;m--;){if(Y(c[m]-a)<=b)return c[m]}else{c=+c;m=a%c;if(m<b)return a-m;if(m>c-b)return a-m+c}return a};a.getRGB=A(function(c){if(!c||(c=J(c)).indexOf("-")+1)return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka};if("none"==c)return{r:-1,g:-1,b:-1,hex:"none",toString:ka};!X[h](c.toLowerCase().substring(0,
2))&&"#"!=c.charAt()&&(c=T(c));if(!c)return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka};var b,m,e,f,d;if(c=c.match(F)){c[2]&&(e=U(c[2].substring(5),16),m=U(c[2].substring(3,5),16),b=U(c[2].substring(1,3),16));c[3]&&(e=U((d=c[3].charAt(3))+d,16),m=U((d=c[3].charAt(2))+d,16),b=U((d=c[3].charAt(1))+d,16));c[4]&&(d=c[4].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b*=2.55),m=K(d[1]),"%"==d[1].slice(-1)&&(m*=2.55),e=K(d[2]),"%"==d[2].slice(-1)&&(e*=2.55),"rgba"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),
d[3]&&"%"==d[3].slice(-1)&&(f/=100));if(c[5])return d=c[5].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b/=100),m=K(d[1]),"%"==d[1].slice(-1)&&(m/=100),e=K(d[2]),"%"==d[2].slice(-1)&&(e/=100),"deg"!=d[0].slice(-3)&&"\u00b0"!=d[0].slice(-1)||(b/=360),"hsba"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),d[3]&&"%"==d[3].slice(-1)&&(f/=100),a.hsb2rgb(b,m,e,f);if(c[6])return d=c[6].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b/=100),m=K(d[1]),"%"==d[1].slice(-1)&&(m/=100),e=K(d[2]),"%"==d[2].slice(-1)&&(e/=100),
"deg"!=d[0].slice(-3)&&"\u00b0"!=d[0].slice(-1)||(b/=360),"hsla"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),d[3]&&"%"==d[3].slice(-1)&&(f/=100),a.hsl2rgb(b,m,e,f);b=Q(I.round(b),255);m=Q(I.round(m),255);e=Q(I.round(e),255);f=Q(P(f,0),1);c={r:b,g:m,b:e,toString:ka};c.hex="#"+(16777216|e|m<<8|b<<16).toString(16).slice(1);c.opacity=y(f,"finite")?f:1;return c}return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka}},a);a.hsb=A(function(c,b,m){return a.hsb2rgb(c,b,m).hex});a.hsl=A(function(c,b,m){return a.hsl2rgb(c,
b,m).hex});a.rgb=A(function(c,a,b,m){if(y(m,"finite")){var e=I.round;return"rgba("+[e(c),e(a),e(b),+m.toFixed(2)]+")"}return"#"+(16777216|b|a<<8|c<<16).toString(16).slice(1)});var T=function(c){var a=G.doc.getElementsByTagName("head")[0]||G.doc.getElementsByTagName("svg")[0];T=A(function(c){if("red"==c.toLowerCase())return"rgb(255, 0, 0)";a.style.color="rgb(255, 0, 0)";a.style.color=c;c=G.doc.defaultView.getComputedStyle(a,aa).getPropertyValue("color");return"rgb(255, 0, 0)"==c?null:c});return T(c)},
qa=function(){return"hsb("+[this.h,this.s,this.b]+")"},ra=function(){return"hsl("+[this.h,this.s,this.l]+")"},ka=function(){return 1==this.opacity||null==this.opacity?this.hex:"rgba("+[this.r,this.g,this.b,this.opacity]+")"},D=function(c,b,m){null==b&&y(c,"object")&&"r"in c&&"g"in c&&"b"in c&&(m=c.b,b=c.g,c=c.r);null==b&&y(c,string)&&(m=a.getRGB(c),c=m.r,b=m.g,m=m.b);if(1<c||1<b||1<m)c/=255,b/=255,m/=255;return[c,b,m]},oa=function(c,b,m,e){c=I.round(255*c);b=I.round(255*b);m=I.round(255*m);c={r:c,
g:b,b:m,opacity:y(e,"finite")?e:1,hex:a.rgb(c,b,m),toString:ka};y(e,"finite")&&(c.opacity=e);return c};a.color=function(c){var b;y(c,"object")&&"h"in c&&"s"in c&&"b"in c?(b=a.hsb2rgb(c),c.r=b.r,c.g=b.g,c.b=b.b,c.opacity=1,c.hex=b.hex):y(c,"object")&&"h"in c&&"s"in c&&"l"in c?(b=a.hsl2rgb(c),c.r=b.r,c.g=b.g,c.b=b.b,c.opacity=1,c.hex=b.hex):(y(c,"string")&&(c=a.getRGB(c)),y(c,"object")&&"r"in c&&"g"in c&&"b"in c&&!("error"in c)?(b=a.rgb2hsl(c),c.h=b.h,c.s=b.s,c.l=b.l,b=a.rgb2hsb(c),c.v=b.b):(c={hex:"none"},
c.r=c.g=c.b=c.h=c.s=c.v=c.l=-1,c.error=1));c.toString=ka;return c};a.hsb2rgb=function(c,a,b,m){y(c,"object")&&"h"in c&&"s"in c&&"b"in c&&(b=c.b,a=c.s,c=c.h,m=c.o);var e,h,d;c=360*c%360/60;d=b*a;a=d*(1-Y(c%2-1));b=e=h=b-d;c=~~c;b+=[d,a,0,0,a,d][c];e+=[a,d,d,a,0,0][c];h+=[0,0,a,d,d,a][c];return oa(b,e,h,m)};a.hsl2rgb=function(c,a,b,m){y(c,"object")&&"h"in c&&"s"in c&&"l"in c&&(b=c.l,a=c.s,c=c.h);if(1<c||1<a||1<b)c/=360,a/=100,b/=100;var e,h,d;c=360*c%360/60;d=2*a*(0.5>b?b:1-b);a=d*(1-Y(c%2-1));b=e=
h=b-d/2;c=~~c;b+=[d,a,0,0,a,d][c];e+=[a,d,d,a,0,0][c];h+=[0,0,a,d,d,a][c];return oa(b,e,h,m)};a.rgb2hsb=function(c,a,b){b=D(c,a,b);c=b[0];a=b[1];b=b[2];var m,e;m=P(c,a,b);e=m-Q(c,a,b);c=((0==e?0:m==c?(a-b)/e:m==a?(b-c)/e+2:(c-a)/e+4)+360)%6*60/360;return{h:c,s:0==e?0:e/m,b:m,toString:qa}};a.rgb2hsl=function(c,a,b){b=D(c,a,b);c=b[0];a=b[1];b=b[2];var m,e,h;m=P(c,a,b);e=Q(c,a,b);h=m-e;c=((0==h?0:m==c?(a-b)/h:m==a?(b-c)/h+2:(c-a)/h+4)+360)%6*60/360;m=(m+e)/2;return{h:c,s:0==h?0:0.5>m?h/(2*m):h/(2-2*
m),l:m,toString:ra}};a.parsePathString=function(c){if(!c)return null;var b=a.path(c);if(b.arr)return a.path.clone(b.arr);var m={a:7,c:6,o:2,h:1,l:2,m:2,r:4,q:4,s:4,t:2,v:1,u:3,z:0},e=[];y(c,"array")&&y(c[0],"array")&&(e=a.path.clone(c));e.length||J(c).replace(W,function(c,a,b){var h=[];c=a.toLowerCase();b.replace(Z,function(c,a){a&&h.push(+a)});"m"==c&&2<h.length&&(e.push([a].concat(h.splice(0,2))),c="l",a="m"==a?"l":"L");"o"==c&&1==h.length&&e.push([a,h[0] ]);if("r"==c)e.push([a].concat(h));else for(;h.length>=
m[c]&&(e.push([a].concat(h.splice(0,m[c]))),m[c]););});e.toString=a.path.toString;b.arr=a.path.clone(e);return e};var O=a.parseTransformString=function(c){if(!c)return null;var b=[];y(c,"array")&&y(c[0],"array")&&(b=a.path.clone(c));b.length||J(c).replace(ma,function(c,a,m){var e=[];a.toLowerCase();m.replace(Z,function(c,a){a&&e.push(+a)});b.push([a].concat(e))});b.toString=a.path.toString;return b};a._.svgTransform2string=d;a._.rgTransform=RegExp("^[a-z][\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*-?\\.?\\d",
"i");a._.transform2matrix=f;a._unit2px=b;a._.getSomeDefs=u;a._.getSomeSVG=p;a.select=function(c){return x(G.doc.querySelector(c))};a.selectAll=function(c){c=G.doc.querySelectorAll(c);for(var b=(a.set||Array)(),m=0;m<c.length;m++)b.push(x(c[m]));return b};setInterval(function(){for(var c in E)if(E[h](c)){var a=E[c],b=a.node;("svg"!=a.type&&!b.ownerSVGElement||"svg"==a.type&&(!b.parentNode||"ownerSVGElement"in b.parentNode&&!b.ownerSVGElement))&&delete E[c]}},1E4);(function(c){function m(c){function a(c,
b){var m=v(c.node,b);(m=(m=m&&m.match(d))&&m[2])&&"#"==m.charAt()&&(m=m.substring(1))&&(f[m]=(f[m]||[]).concat(function(a){var m={};m[b]=ca(a);v(c.node,m)}))}function b(c){var a=v(c.node,"xlink:href");a&&"#"==a.charAt()&&(a=a.substring(1))&&(f[a]=(f[a]||[]).concat(function(a){c.attr("xlink:href","#"+a)}))}var e=c.selectAll("*"),h,d=/^\s*url\(("|'|)(.*)\1\)\s*$/;c=[];for(var f={},l=0,E=e.length;l<E;l++){h=e[l];a(h,"fill");a(h,"stroke");a(h,"filter");a(h,"mask");a(h,"clip-path");b(h);var t=v(h.node,
"id");t&&(v(h.node,{id:h.id}),c.push({old:t,id:h.id}))}l=0;for(E=c.length;l<E;l++)if(e=f[c[l].old])for(h=0,t=e.length;h<t;h++)e[h](c[l].id)}function e(c,a,b){return function(m){m=m.slice(c,a);1==m.length&&(m=m[0]);return b?b(m):m}}function d(c){return function(){var a=c?"<"+this.type:"",b=this.node.attributes,m=this.node.childNodes;if(c)for(var e=0,h=b.length;e<h;e++)a+=" "+b[e].name+'="'+b[e].value.replace(/"/g,'\\"')+'"';if(m.length){c&&(a+=">");e=0;for(h=m.length;e<h;e++)3==m[e].nodeType?a+=m[e].nodeValue:
1==m[e].nodeType&&(a+=x(m[e]).toString());c&&(a+="</"+this.type+">")}else c&&(a+="/>");return a}}c.attr=function(c,a){if(!c)return this;if(y(c,"string"))if(1<arguments.length){var b={};b[c]=a;c=b}else return k("snap.util.getattr."+c,this).firstDefined();for(var m in c)c[h](m)&&k("snap.util.attr."+m,this,c[m]);return this};c.getBBox=function(c){if(!a.Matrix||!a.path)return this.node.getBBox();var b=this,m=new a.Matrix;if(b.removed)return a._.box();for(;"use"==b.type;)if(c||(m=m.add(b.transform().localMatrix.translate(b.attr("x")||
0,b.attr("y")||0))),b.original)b=b.original;else var e=b.attr("xlink:href"),b=b.original=b.node.ownerDocument.getElementById(e.substring(e.indexOf("#")+1));var e=b._,h=a.path.get[b.type]||a.path.get.deflt;try{if(c)return e.bboxwt=h?a.path.getBBox(b.realPath=h(b)):a._.box(b.node.getBBox()),a._.box(e.bboxwt);b.realPath=h(b);b.matrix=b.transform().localMatrix;e.bbox=a.path.getBBox(a.path.map(b.realPath,m.add(b.matrix)));return a._.box(e.bbox)}catch(d){return a._.box()}};var f=function(){return this.string};
c.transform=function(c){var b=this._;if(null==c){var m=this;c=new a.Matrix(this.node.getCTM());for(var e=n(this),h=[e],d=new a.Matrix,l=e.toTransformString(),b=J(e)==J(this.matrix)?J(b.transform):l;"svg"!=m.type&&(m=m.parent());)h.push(n(m));for(m=h.length;m--;)d.add(h[m]);return{string:b,globalMatrix:c,totalMatrix:d,localMatrix:e,diffMatrix:c.clone().add(e.invert()),global:c.toTransformString(),total:d.toTransformString(),local:l,toString:f}}c instanceof a.Matrix?this.matrix=c:n(this,c);this.node&&
("linearGradient"==this.type||"radialGradient"==this.type?v(this.node,{gradientTransform:this.matrix}):"pattern"==this.type?v(this.node,{patternTransform:this.matrix}):v(this.node,{transform:this.matrix}));return this};c.parent=function(){return x(this.node.parentNode)};c.append=c.add=function(c){if(c){if("set"==c.type){var a=this;c.forEach(function(c){a.add(c)});return this}c=x(c);this.node.appendChild(c.node);c.paper=this.paper}return this};c.appendTo=function(c){c&&(c=x(c),c.append(this));return this};
c.prepend=function(c){if(c){if("set"==c.type){var a=this,b;c.forEach(function(c){b?b.after(c):a.prepend(c);b=c});return this}c=x(c);var m=c.parent();this.node.insertBefore(c.node,this.node.firstChild);this.add&&this.add();c.paper=this.paper;this.parent()&&this.parent().add();m&&m.add()}return this};c.prependTo=function(c){c=x(c);c.prepend(this);return this};c.before=function(c){if("set"==c.type){var a=this;c.forEach(function(c){var b=c.parent();a.node.parentNode.insertBefore(c.node,a.node);b&&b.add()});
this.parent().add();return this}c=x(c);var b=c.parent();this.node.parentNode.insertBefore(c.node,this.node);this.parent()&&this.parent().add();b&&b.add();c.paper=this.paper;return this};c.after=function(c){c=x(c);var a=c.parent();this.node.nextSibling?this.node.parentNode.insertBefore(c.node,this.node.nextSibling):this.node.parentNode.appendChild(c.node);this.parent()&&this.parent().add();a&&a.add();c.paper=this.paper;return this};c.insertBefore=function(c){c=x(c);var a=this.parent();c.node.parentNode.insertBefore(this.node,
c.node);this.paper=c.paper;a&&a.add();c.parent()&&c.parent().add();return this};c.insertAfter=function(c){c=x(c);var a=this.parent();c.node.parentNode.insertBefore(this.node,c.node.nextSibling);this.paper=c.paper;a&&a.add();c.parent()&&c.parent().add();return this};c.remove=function(){var c=this.parent();this.node.parentNode&&this.node.parentNode.removeChild(this.node);delete this.paper;this.removed=!0;c&&c.add();return this};c.select=function(c){return x(this.node.querySelector(c))};c.selectAll=
function(c){c=this.node.querySelectorAll(c);for(var b=(a.set||Array)(),m=0;m<c.length;m++)b.push(x(c[m]));return b};c.asPX=function(c,a){null==a&&(a=this.attr(c));return+b(this,c,a)};c.use=function(){var c,a=this.node.id;a||(a=this.id,v(this.node,{id:a}));c="linearGradient"==this.type||"radialGradient"==this.type||"pattern"==this.type?r(this.type,this.node.parentNode):r("use",this.node.parentNode);v(c.node,{"xlink:href":"#"+a});c.original=this;return c};var l=/\S+/g;c.addClass=function(c){var a=(c||
"").match(l)||[];c=this.node;var b=c.className.baseVal,m=b.match(l)||[],e,h,d;if(a.length){for(e=0;d=a[e++];)h=m.indexOf(d),~h||m.push(d);a=m.join(" ");b!=a&&(c.className.baseVal=a)}return this};c.removeClass=function(c){var a=(c||"").match(l)||[];c=this.node;var b=c.className.baseVal,m=b.match(l)||[],e,h;if(m.length){for(e=0;h=a[e++];)h=m.indexOf(h),~h&&m.splice(h,1);a=m.join(" ");b!=a&&(c.className.baseVal=a)}return this};c.hasClass=function(c){return!!~(this.node.className.baseVal.match(l)||[]).indexOf(c)};
c.toggleClass=function(c,a){if(null!=a)return a?this.addClass(c):this.removeClass(c);var b=(c||"").match(l)||[],m=this.node,e=m.className.baseVal,h=e.match(l)||[],d,f,E;for(d=0;E=b[d++];)f=h.indexOf(E),~f?h.splice(f,1):h.push(E);b=h.join(" ");e!=b&&(m.className.baseVal=b);return this};c.clone=function(){var c=x(this.node.cloneNode(!0));v(c.node,"id")&&v(c.node,{id:c.id});m(c);c.insertAfter(this);return c};c.toDefs=function(){u(this).appendChild(this.node);return this};c.pattern=c.toPattern=function(c,
a,b,m){var e=r("pattern",u(this));null==c&&(c=this.getBBox());y(c,"object")&&"x"in c&&(a=c.y,b=c.width,m=c.height,c=c.x);v(e.node,{x:c,y:a,width:b,height:m,patternUnits:"userSpaceOnUse",id:e.id,viewBox:[c,a,b,m].join(" ")});e.node.appendChild(this.node);return e};c.marker=function(c,a,b,m,e,h){var d=r("marker",u(this));null==c&&(c=this.getBBox());y(c,"object")&&"x"in c&&(a=c.y,b=c.width,m=c.height,e=c.refX||c.cx,h=c.refY||c.cy,c=c.x);v(d.node,{viewBox:[c,a,b,m].join(" "),markerWidth:b,markerHeight:m,
orient:"auto",refX:e||0,refY:h||0,id:d.id});d.node.appendChild(this.node);return d};var E=function(c,a,b,m){"function"!=typeof b||b.length||(m=b,b=L.linear);this.attr=c;this.dur=a;b&&(this.easing=b);m&&(this.callback=m)};a._.Animation=E;a.animation=function(c,a,b,m){return new E(c,a,b,m)};c.inAnim=function(){var c=[],a;for(a in this.anims)this.anims[h](a)&&function(a){c.push({anim:new E(a._attrs,a.dur,a.easing,a._callback),mina:a,curStatus:a.status(),status:function(c){return a.status(c)},stop:function(){a.stop()}})}(this.anims[a]);
return c};a.animate=function(c,a,b,m,e,h){"function"!=typeof e||e.length||(h=e,e=L.linear);var d=L.time();c=L(c,a,d,d+m,L.time,b,e);h&&k.once("mina.finish."+c.id,h);return c};c.stop=function(){for(var c=this.inAnim(),a=0,b=c.length;a<b;a++)c[a].stop();return this};c.animate=function(c,a,b,m){"function"!=typeof b||b.length||(m=b,b=L.linear);c instanceof E&&(m=c.callback,b=c.easing,a=b.dur,c=c.attr);var d=[],f=[],l={},t,ca,n,T=this,q;for(q in c)if(c[h](q)){T.equal?(n=T.equal(q,J(c[q])),t=n.from,ca=
n.to,n=n.f):(t=+T.attr(q),ca=+c[q]);var la=y(t,"array")?t.length:1;l[q]=e(d.length,d.length+la,n);d=d.concat(t);f=f.concat(ca)}t=L.time();var p=L(d,f,t,t+a,L.time,function(c){var a={},b;for(b in l)l[h](b)&&(a[b]=l[b](c));T.attr(a)},b);T.anims[p.id]=p;p._attrs=c;p._callback=m;k("snap.animcreated."+T.id,p);k.once("mina.finish."+p.id,function(){delete T.anims[p.id];m&&m.call(T)});k.once("mina.stop."+p.id,function(){delete T.anims[p.id]});return T};var T={};c.data=function(c,b){var m=T[this.id]=T[this.id]||
{};if(0==arguments.length)return k("snap.data.get."+this.id,this,m,null),m;if(1==arguments.length){if(a.is(c,"object")){for(var e in c)c[h](e)&&this.data(e,c[e]);return this}k("snap.data.get."+this.id,this,m[c],c);return m[c]}m[c]=b;k("snap.data.set."+this.id,this,b,c);return this};c.removeData=function(c){null==c?T[this.id]={}:T[this.id]&&delete T[this.id][c];return this};c.outerSVG=c.toString=d(1);c.innerSVG=d()})(e.prototype);a.parse=function(c){var a=G.doc.createDocumentFragment(),b=!0,m=G.doc.createElement("div");
c=J(c);c.match(/^\s*<\s*svg(?:\s|>)/)||(c="<svg>"+c+"</svg>",b=!1);m.innerHTML=c;if(c=m.getElementsByTagName("svg")[0])if(b)a=c;else for(;c.firstChild;)a.appendChild(c.firstChild);m.innerHTML=aa;return new l(a)};l.prototype.select=e.prototype.select;l.prototype.selectAll=e.prototype.selectAll;a.fragment=function(){for(var c=Array.prototype.slice.call(arguments,0),b=G.doc.createDocumentFragment(),m=0,e=c.length;m<e;m++){var h=c[m];h.node&&h.node.nodeType&&b.appendChild(h.node);h.nodeType&&b.appendChild(h);
"string"==typeof h&&b.appendChild(a.parse(h).node)}return new l(b)};a._.make=r;a._.wrap=x;s.prototype.el=function(c,a){var b=r(c,this.node);a&&b.attr(a);return b};k.on("snap.util.getattr",function(){var c=k.nt(),c=c.substring(c.lastIndexOf(".")+1),a=c.replace(/[A-Z]/g,function(c){return"-"+c.toLowerCase()});return pa[h](a)?this.node.ownerDocument.defaultView.getComputedStyle(this.node,null).getPropertyValue(a):v(this.node,c)});var pa={"alignment-baseline":0,"baseline-shift":0,clip:0,"clip-path":0,
"clip-rule":0,color:0,"color-interpolation":0,"color-interpolation-filters":0,"color-profile":0,"color-rendering":0,cursor:0,direction:0,display:0,"dominant-baseline":0,"enable-background":0,fill:0,"fill-opacity":0,"fill-rule":0,filter:0,"flood-color":0,"flood-opacity":0,font:0,"font-family":0,"font-size":0,"font-size-adjust":0,"font-stretch":0,"font-style":0,"font-variant":0,"font-weight":0,"glyph-orientation-horizontal":0,"glyph-orientation-vertical":0,"image-rendering":0,kerning:0,"letter-spacing":0,
"lighting-color":0,marker:0,"marker-end":0,"marker-mid":0,"marker-start":0,mask:0,opacity:0,overflow:0,"pointer-events":0,"shape-rendering":0,"stop-color":0,"stop-opacity":0,stroke:0,"stroke-dasharray":0,"stroke-dashoffset":0,"stroke-linecap":0,"stroke-linejoin":0,"stroke-miterlimit":0,"stroke-opacity":0,"stroke-width":0,"text-anchor":0,"text-decoration":0,"text-rendering":0,"unicode-bidi":0,visibility:0,"word-spacing":0,"writing-mode":0};k.on("snap.util.attr",function(c){var a=k.nt(),b={},a=a.substring(a.lastIndexOf(".")+
1);b[a]=c;var m=a.replace(/-(\w)/gi,function(c,a){return a.toUpperCase()}),a=a.replace(/[A-Z]/g,function(c){return"-"+c.toLowerCase()});pa[h](a)?this.node.style[m]=null==c?aa:c:v(this.node,b)});a.ajax=function(c,a,b,m){var e=new XMLHttpRequest,h=V();if(e){if(y(a,"function"))m=b,b=a,a=null;else if(y(a,"object")){var d=[],f;for(f in a)a.hasOwnProperty(f)&&d.push(encodeURIComponent(f)+"="+encodeURIComponent(a[f]));a=d.join("&")}e.open(a?"POST":"GET",c,!0);a&&(e.setRequestHeader("X-Requested-With","XMLHttpRequest"),
e.setRequestHeader("Content-type","application/x-www-form-urlencoded"));b&&(k.once("snap.ajax."+h+".0",b),k.once("snap.ajax."+h+".200",b),k.once("snap.ajax."+h+".304",b));e.onreadystatechange=function(){4==e.readyState&&k("snap.ajax."+h+"."+e.status,m,e)};if(4==e.readyState)return e;e.send(a);return e}};a.load=function(c,b,m){a.ajax(c,function(c){c=a.parse(c.responseText);m?b.call(m,c):b(c)})};a.getElementByPoint=function(c,a){var b,m,e=G.doc.elementFromPoint(c,a);if(G.win.opera&&"svg"==e.tagName){b=
e;m=b.getBoundingClientRect();b=b.ownerDocument;var h=b.body,d=b.documentElement;b=m.top+(g.win.pageYOffset||d.scrollTop||h.scrollTop)-(d.clientTop||h.clientTop||0);m=m.left+(g.win.pageXOffset||d.scrollLeft||h.scrollLeft)-(d.clientLeft||h.clientLeft||0);h=e.createSVGRect();h.x=c-m;h.y=a-b;h.width=h.height=1;b=e.getIntersectionList(h,null);b.length&&(e=b[b.length-1])}return e?x(e):null};a.plugin=function(c){c(a,e,s,G,l)};return G.win.Snap=a}();C.plugin(function(a,k,y,M,A){function w(a,d,f,b,q,e){null==
d&&"[object SVGMatrix]"==z.call(a)?(this.a=a.a,this.b=a.b,this.c=a.c,this.d=a.d,this.e=a.e,this.f=a.f):null!=a?(this.a=+a,this.b=+d,this.c=+f,this.d=+b,this.e=+q,this.f=+e):(this.a=1,this.c=this.b=0,this.d=1,this.f=this.e=0)}var z=Object.prototype.toString,d=String,f=Math;(function(n){function k(a){return a[0]*a[0]+a[1]*a[1]}function p(a){var d=f.sqrt(k(a));a[0]&&(a[0]/=d);a[1]&&(a[1]/=d)}n.add=function(a,d,e,f,n,p){var k=[[],[],[] ],u=[[this.a,this.c,this.e],[this.b,this.d,this.f],[0,0,1] ];d=[[a,
e,n],[d,f,p],[0,0,1] ];a&&a instanceof w&&(d=[[a.a,a.c,a.e],[a.b,a.d,a.f],[0,0,1] ]);for(a=0;3>a;a++)for(e=0;3>e;e++){for(f=n=0;3>f;f++)n+=u[a][f]*d[f][e];k[a][e]=n}this.a=k[0][0];this.b=k[1][0];this.c=k[0][1];this.d=k[1][1];this.e=k[0][2];this.f=k[1][2];return this};n.invert=function(){var a=this.a*this.d-this.b*this.c;return new w(this.d/a,-this.b/a,-this.c/a,this.a/a,(this.c*this.f-this.d*this.e)/a,(this.b*this.e-this.a*this.f)/a)};n.clone=function(){return new w(this.a,this.b,this.c,this.d,this.e,
this.f)};n.translate=function(a,d){return this.add(1,0,0,1,a,d)};n.scale=function(a,d,e,f){null==d&&(d=a);(e||f)&&this.add(1,0,0,1,e,f);this.add(a,0,0,d,0,0);(e||f)&&this.add(1,0,0,1,-e,-f);return this};n.rotate=function(b,d,e){b=a.rad(b);d=d||0;e=e||0;var l=+f.cos(b).toFixed(9);b=+f.sin(b).toFixed(9);this.add(l,b,-b,l,d,e);return this.add(1,0,0,1,-d,-e)};n.x=function(a,d){return a*this.a+d*this.c+this.e};n.y=function(a,d){return a*this.b+d*this.d+this.f};n.get=function(a){return+this[d.fromCharCode(97+
a)].toFixed(4)};n.toString=function(){return"matrix("+[this.get(0),this.get(1),this.get(2),this.get(3),this.get(4),this.get(5)].join()+")"};n.offset=function(){return[this.e.toFixed(4),this.f.toFixed(4)]};n.determinant=function(){return this.a*this.d-this.b*this.c};n.split=function(){var b={};b.dx=this.e;b.dy=this.f;var d=[[this.a,this.c],[this.b,this.d] ];b.scalex=f.sqrt(k(d[0]));p(d[0]);b.shear=d[0][0]*d[1][0]+d[0][1]*d[1][1];d[1]=[d[1][0]-d[0][0]*b.shear,d[1][1]-d[0][1]*b.shear];b.scaley=f.sqrt(k(d[1]));
p(d[1]);b.shear/=b.scaley;0>this.determinant()&&(b.scalex=-b.scalex);var e=-d[0][1],d=d[1][1];0>d?(b.rotate=a.deg(f.acos(d)),0>e&&(b.rotate=360-b.rotate)):b.rotate=a.deg(f.asin(e));b.isSimple=!+b.shear.toFixed(9)&&(b.scalex.toFixed(9)==b.scaley.toFixed(9)||!b.rotate);b.isSuperSimple=!+b.shear.toFixed(9)&&b.scalex.toFixed(9)==b.scaley.toFixed(9)&&!b.rotate;b.noRotation=!+b.shear.toFixed(9)&&!b.rotate;return b};n.toTransformString=function(a){a=a||this.split();if(+a.shear.toFixed(9))return"m"+[this.get(0),
this.get(1),this.get(2),this.get(3),this.get(4),this.get(5)];a.scalex=+a.scalex.toFixed(4);a.scaley=+a.scaley.toFixed(4);a.rotate=+a.rotate.toFixed(4);return(a.dx||a.dy?"t"+[+a.dx.toFixed(4),+a.dy.toFixed(4)]:"")+(1!=a.scalex||1!=a.scaley?"s"+[a.scalex,a.scaley,0,0]:"")+(a.rotate?"r"+[+a.rotate.toFixed(4),0,0]:"")}})(w.prototype);a.Matrix=w;a.matrix=function(a,d,f,b,k,e){return new w(a,d,f,b,k,e)}});C.plugin(function(a,v,y,M,A){function w(h){return function(d){k.stop();d instanceof A&&1==d.node.childNodes.length&&
("radialGradient"==d.node.firstChild.tagName||"linearGradient"==d.node.firstChild.tagName||"pattern"==d.node.firstChild.tagName)&&(d=d.node.firstChild,b(this).appendChild(d),d=u(d));if(d instanceof v)if("radialGradient"==d.type||"linearGradient"==d.type||"pattern"==d.type){d.node.id||e(d.node,{id:d.id});var f=l(d.node.id)}else f=d.attr(h);else f=a.color(d),f.error?(f=a(b(this).ownerSVGElement).gradient(d))?(f.node.id||e(f.node,{id:f.id}),f=l(f.node.id)):f=d:f=r(f);d={};d[h]=f;e(this.node,d);this.node.style[h]=
x}}function z(a){k.stop();a==+a&&(a+="px");this.node.style.fontSize=a}function d(a){var b=[];a=a.childNodes;for(var e=0,f=a.length;e<f;e++){var l=a[e];3==l.nodeType&&b.push(l.nodeValue);"tspan"==l.tagName&&(1==l.childNodes.length&&3==l.firstChild.nodeType?b.push(l.firstChild.nodeValue):b.push(d(l)))}return b}function f(){k.stop();return this.node.style.fontSize}var n=a._.make,u=a._.wrap,p=a.is,b=a._.getSomeDefs,q=/^url\(#?([^)]+)\)$/,e=a._.$,l=a.url,r=String,s=a._.separator,x="";k.on("snap.util.attr.mask",
function(a){if(a instanceof v||a instanceof A){k.stop();a instanceof A&&1==a.node.childNodes.length&&(a=a.node.firstChild,b(this).appendChild(a),a=u(a));if("mask"==a.type)var d=a;else d=n("mask",b(this)),d.node.appendChild(a.node);!d.node.id&&e(d.node,{id:d.id});e(this.node,{mask:l(d.id)})}});(function(a){k.on("snap.util.attr.clip",a);k.on("snap.util.attr.clip-path",a);k.on("snap.util.attr.clipPath",a)})(function(a){if(a instanceof v||a instanceof A){k.stop();if("clipPath"==a.type)var d=a;else d=
n("clipPath",b(this)),d.node.appendChild(a.node),!d.node.id&&e(d.node,{id:d.id});e(this.node,{"clip-path":l(d.id)})}});k.on("snap.util.attr.fill",w("fill"));k.on("snap.util.attr.stroke",w("stroke"));var G=/^([lr])(?:\(([^)]*)\))?(.*)$/i;k.on("snap.util.grad.parse",function(a){a=r(a);var b=a.match(G);if(!b)return null;a=b[1];var e=b[2],b=b[3],e=e.split(/\s*,\s*/).map(function(a){return+a==a?+a:a});1==e.length&&0==e[0]&&(e=[]);b=b.split("-");b=b.map(function(a){a=a.split(":");var b={color:a[0]};a[1]&&
(b.offset=parseFloat(a[1]));return b});return{type:a,params:e,stops:b}});k.on("snap.util.attr.d",function(b){k.stop();p(b,"array")&&p(b[0],"array")&&(b=a.path.toString.call(b));b=r(b);b.match(/[ruo]/i)&&(b=a.path.toAbsolute(b));e(this.node,{d:b})})(-1);k.on("snap.util.attr.#text",function(a){k.stop();a=r(a);for(a=M.doc.createTextNode(a);this.node.firstChild;)this.node.removeChild(this.node.firstChild);this.node.appendChild(a)})(-1);k.on("snap.util.attr.path",function(a){k.stop();this.attr({d:a})})(-1);
k.on("snap.util.attr.class",function(a){k.stop();this.node.className.baseVal=a})(-1);k.on("snap.util.attr.viewBox",function(a){a=p(a,"object")&&"x"in a?[a.x,a.y,a.width,a.height].join(" "):p(a,"array")?a.join(" "):a;e(this.node,{viewBox:a});k.stop()})(-1);k.on("snap.util.attr.transform",function(a){this.transform(a);k.stop()})(-1);k.on("snap.util.attr.r",function(a){"rect"==this.type&&(k.stop(),e(this.node,{rx:a,ry:a}))})(-1);k.on("snap.util.attr.textpath",function(a){k.stop();if("text"==this.type){var d,
f;if(!a&&this.textPath){for(a=this.textPath;a.node.firstChild;)this.node.appendChild(a.node.firstChild);a.remove();delete this.textPath}else if(p(a,"string")?(d=b(this),a=u(d.parentNode).path(a),d.appendChild(a.node),d=a.id,a.attr({id:d})):(a=u(a),a instanceof v&&(d=a.attr("id"),d||(d=a.id,a.attr({id:d})))),d)if(a=this.textPath,f=this.node,a)a.attr({"xlink:href":"#"+d});else{for(a=e("textPath",{"xlink:href":"#"+d});f.firstChild;)a.appendChild(f.firstChild);f.appendChild(a);this.textPath=u(a)}}})(-1);
k.on("snap.util.attr.text",function(a){if("text"==this.type){for(var b=this.node,d=function(a){var b=e("tspan");if(p(a,"array"))for(var f=0;f<a.length;f++)b.appendChild(d(a[f]));else b.appendChild(M.doc.createTextNode(a));b.normalize&&b.normalize();return b};b.firstChild;)b.removeChild(b.firstChild);for(a=d(a);a.firstChild;)b.appendChild(a.firstChild)}k.stop()})(-1);k.on("snap.util.attr.fontSize",z)(-1);k.on("snap.util.attr.font-size",z)(-1);k.on("snap.util.getattr.transform",function(){k.stop();
return this.transform()})(-1);k.on("snap.util.getattr.textpath",function(){k.stop();return this.textPath})(-1);(function(){function b(d){return function(){k.stop();var b=M.doc.defaultView.getComputedStyle(this.node,null).getPropertyValue("marker-"+d);return"none"==b?b:a(M.doc.getElementById(b.match(q)[1]))}}function d(a){return function(b){k.stop();var d="marker"+a.charAt(0).toUpperCase()+a.substring(1);if(""==b||!b)this.node.style[d]="none";else if("marker"==b.type){var f=b.node.id;f||e(b.node,{id:b.id});
this.node.style[d]=l(f)}}}k.on("snap.util.getattr.marker-end",b("end"))(-1);k.on("snap.util.getattr.markerEnd",b("end"))(-1);k.on("snap.util.getattr.marker-start",b("start"))(-1);k.on("snap.util.getattr.markerStart",b("start"))(-1);k.on("snap.util.getattr.marker-mid",b("mid"))(-1);k.on("snap.util.getattr.markerMid",b("mid"))(-1);k.on("snap.util.attr.marker-end",d("end"))(-1);k.on("snap.util.attr.markerEnd",d("end"))(-1);k.on("snap.util.attr.marker-start",d("start"))(-1);k.on("snap.util.attr.markerStart",
d("start"))(-1);k.on("snap.util.attr.marker-mid",d("mid"))(-1);k.on("snap.util.attr.markerMid",d("mid"))(-1)})();k.on("snap.util.getattr.r",function(){if("rect"==this.type&&e(this.node,"rx")==e(this.node,"ry"))return k.stop(),e(this.node,"rx")})(-1);k.on("snap.util.getattr.text",function(){if("text"==this.type||"tspan"==this.type){k.stop();var a=d(this.node);return 1==a.length?a[0]:a}})(-1);k.on("snap.util.getattr.#text",function(){return this.node.textContent})(-1);k.on("snap.util.getattr.viewBox",
function(){k.stop();var b=e(this.node,"viewBox");if(b)return b=b.split(s),a._.box(+b[0],+b[1],+b[2],+b[3])})(-1);k.on("snap.util.getattr.points",function(){var a=e(this.node,"points");k.stop();if(a)return a.split(s)})(-1);k.on("snap.util.getattr.path",function(){var a=e(this.node,"d");k.stop();return a})(-1);k.on("snap.util.getattr.class",function(){return this.node.className.baseVal})(-1);k.on("snap.util.getattr.fontSize",f)(-1);k.on("snap.util.getattr.font-size",f)(-1)});C.plugin(function(a,v,y,
M,A){function w(a){return a}function z(a){return function(b){return+b.toFixed(3)+a}}var d={"+":function(a,b){return a+b},"-":function(a,b){return a-b},"/":function(a,b){return a/b},"*":function(a,b){return a*b}},f=String,n=/[a-z]+$/i,u=/^\s*([+\-\/*])\s*=\s*([\d.eE+\-]+)\s*([^\d\s]+)?\s*$/;k.on("snap.util.attr",function(a){if(a=f(a).match(u)){var b=k.nt(),b=b.substring(b.lastIndexOf(".")+1),q=this.attr(b),e={};k.stop();var l=a[3]||"",r=q.match(n),s=d[a[1] ];r&&r==l?a=s(parseFloat(q),+a[2]):(q=this.asPX(b),
a=s(this.asPX(b),this.asPX(b,a[2]+l)));isNaN(q)||isNaN(a)||(e[b]=a,this.attr(e))}})(-10);k.on("snap.util.equal",function(a,b){var q=f(this.attr(a)||""),e=f(b).match(u);if(e){k.stop();var l=e[3]||"",r=q.match(n),s=d[e[1] ];if(r&&r==l)return{from:parseFloat(q),to:s(parseFloat(q),+e[2]),f:z(r)};q=this.asPX(a);return{from:q,to:s(q,this.asPX(a,e[2]+l)),f:w}}})(-10)});C.plugin(function(a,v,y,M,A){var w=y.prototype,z=a.is;w.rect=function(a,d,k,p,b,q){var e;null==q&&(q=b);z(a,"object")&&"[object Object]"==
a?e=a:null!=a&&(e={x:a,y:d,width:k,height:p},null!=b&&(e.rx=b,e.ry=q));return this.el("rect",e)};w.circle=function(a,d,k){var p;z(a,"object")&&"[object Object]"==a?p=a:null!=a&&(p={cx:a,cy:d,r:k});return this.el("circle",p)};var d=function(){function a(){this.parentNode.removeChild(this)}return function(d,k){var p=M.doc.createElement("img"),b=M.doc.body;p.style.cssText="position:absolute;left:-9999em;top:-9999em";p.onload=function(){k.call(p);p.onload=p.onerror=null;b.removeChild(p)};p.onerror=a;
b.appendChild(p);p.src=d}}();w.image=function(f,n,k,p,b){var q=this.el("image");if(z(f,"object")&&"src"in f)q.attr(f);else if(null!=f){var e={"xlink:href":f,preserveAspectRatio:"none"};null!=n&&null!=k&&(e.x=n,e.y=k);null!=p&&null!=b?(e.width=p,e.height=b):d(f,function(){a._.$(q.node,{width:this.offsetWidth,height:this.offsetHeight})});a._.$(q.node,e)}return q};w.ellipse=function(a,d,k,p){var b;z(a,"object")&&"[object Object]"==a?b=a:null!=a&&(b={cx:a,cy:d,rx:k,ry:p});return this.el("ellipse",b)};
w.path=function(a){var d;z(a,"object")&&!z(a,"array")?d=a:a&&(d={d:a});return this.el("path",d)};w.group=w.g=function(a){var d=this.el("g");1==arguments.length&&a&&!a.type?d.attr(a):arguments.length&&d.add(Array.prototype.slice.call(arguments,0));return d};w.svg=function(a,d,k,p,b,q,e,l){var r={};z(a,"object")&&null==d?r=a:(null!=a&&(r.x=a),null!=d&&(r.y=d),null!=k&&(r.width=k),null!=p&&(r.height=p),null!=b&&null!=q&&null!=e&&null!=l&&(r.viewBox=[b,q,e,l]));return this.el("svg",r)};w.mask=function(a){var d=
this.el("mask");1==arguments.length&&a&&!a.type?d.attr(a):arguments.length&&d.add(Array.prototype.slice.call(arguments,0));return d};w.ptrn=function(a,d,k,p,b,q,e,l){if(z(a,"object"))var r=a;else arguments.length?(r={},null!=a&&(r.x=a),null!=d&&(r.y=d),null!=k&&(r.width=k),null!=p&&(r.height=p),null!=b&&null!=q&&null!=e&&null!=l&&(r.viewBox=[b,q,e,l])):r={patternUnits:"userSpaceOnUse"};return this.el("pattern",r)};w.use=function(a){return null!=a?(make("use",this.node),a instanceof v&&(a.attr("id")||
a.attr({id:ID()}),a=a.attr("id")),this.el("use",{"xlink:href":a})):v.prototype.use.call(this)};w.text=function(a,d,k){var p={};z(a,"object")?p=a:null!=a&&(p={x:a,y:d,text:k||""});return this.el("text",p)};w.line=function(a,d,k,p){var b={};z(a,"object")?b=a:null!=a&&(b={x1:a,x2:k,y1:d,y2:p});return this.el("line",b)};w.polyline=function(a){1<arguments.length&&(a=Array.prototype.slice.call(arguments,0));var d={};z(a,"object")&&!z(a,"array")?d=a:null!=a&&(d={points:a});return this.el("polyline",d)};
w.polygon=function(a){1<arguments.length&&(a=Array.prototype.slice.call(arguments,0));var d={};z(a,"object")&&!z(a,"array")?d=a:null!=a&&(d={points:a});return this.el("polygon",d)};(function(){function d(){return this.selectAll("stop")}function n(b,d){var f=e("stop"),k={offset:+d+"%"};b=a.color(b);k["stop-color"]=b.hex;1>b.opacity&&(k["stop-opacity"]=b.opacity);e(f,k);this.node.appendChild(f);return this}function u(){if("linearGradient"==this.type){var b=e(this.node,"x1")||0,d=e(this.node,"x2")||
1,f=e(this.node,"y1")||0,k=e(this.node,"y2")||0;return a._.box(b,f,math.abs(d-b),math.abs(k-f))}b=this.node.r||0;return a._.box((this.node.cx||0.5)-b,(this.node.cy||0.5)-b,2*b,2*b)}function p(a,d){function f(a,b){for(var d=(b-u)/(a-w),e=w;e<a;e++)h[e].offset=+(+u+d*(e-w)).toFixed(2);w=a;u=b}var n=k("snap.util.grad.parse",null,d).firstDefined(),p;if(!n)return null;n.params.unshift(a);p="l"==n.type.toLowerCase()?b.apply(0,n.params):q.apply(0,n.params);n.type!=n.type.toLowerCase()&&e(p.node,{gradientUnits:"userSpaceOnUse"});
var h=n.stops,n=h.length,u=0,w=0;n--;for(var v=0;v<n;v++)"offset"in h[v]&&f(v,h[v].offset);h[n].offset=h[n].offset||100;f(n,h[n].offset);for(v=0;v<=n;v++){var y=h[v];p.addStop(y.color,y.offset)}return p}function b(b,k,p,q,w){b=a._.make("linearGradient",b);b.stops=d;b.addStop=n;b.getBBox=u;null!=k&&e(b.node,{x1:k,y1:p,x2:q,y2:w});return b}function q(b,k,p,q,w,h){b=a._.make("radialGradient",b);b.stops=d;b.addStop=n;b.getBBox=u;null!=k&&e(b.node,{cx:k,cy:p,r:q});null!=w&&null!=h&&e(b.node,{fx:w,fy:h});
return b}var e=a._.$;w.gradient=function(a){return p(this.defs,a)};w.gradientLinear=function(a,d,e,f){return b(this.defs,a,d,e,f)};w.gradientRadial=function(a,b,d,e,f){return q(this.defs,a,b,d,e,f)};w.toString=function(){var b=this.node.ownerDocument,d=b.createDocumentFragment(),b=b.createElement("div"),e=this.node.cloneNode(!0);d.appendChild(b);b.appendChild(e);a._.$(e,{xmlns:"http://www.w3.org/2000/svg"});b=b.innerHTML;d.removeChild(d.firstChild);return b};w.clear=function(){for(var a=this.node.firstChild,
b;a;)b=a.nextSibling,"defs"!=a.tagName?a.parentNode.removeChild(a):w.clear.call({node:a}),a=b}})()});C.plugin(function(a,k,y,M){function A(a){var b=A.ps=A.ps||{};b[a]?b[a].sleep=100:b[a]={sleep:100};setTimeout(function(){for(var d in b)b[L](d)&&d!=a&&(b[d].sleep--,!b[d].sleep&&delete b[d])});return b[a]}function w(a,b,d,e){null==a&&(a=b=d=e=0);null==b&&(b=a.y,d=a.width,e=a.height,a=a.x);return{x:a,y:b,width:d,w:d,height:e,h:e,x2:a+d,y2:b+e,cx:a+d/2,cy:b+e/2,r1:F.min(d,e)/2,r2:F.max(d,e)/2,r0:F.sqrt(d*
d+e*e)/2,path:s(a,b,d,e),vb:[a,b,d,e].join(" ")}}function z(){return this.join(",").replace(N,"$1")}function d(a){a=C(a);a.toString=z;return a}function f(a,b,d,h,f,k,l,n,p){if(null==p)return e(a,b,d,h,f,k,l,n);if(0>p||e(a,b,d,h,f,k,l,n)<p)p=void 0;else{var q=0.5,O=1-q,s;for(s=e(a,b,d,h,f,k,l,n,O);0.01<Z(s-p);)q/=2,O+=(s<p?1:-1)*q,s=e(a,b,d,h,f,k,l,n,O);p=O}return u(a,b,d,h,f,k,l,n,p)}function n(b,d){function e(a){return+(+a).toFixed(3)}return a._.cacher(function(a,h,l){a instanceof k&&(a=a.attr("d"));
a=I(a);for(var n,p,D,q,O="",s={},c=0,t=0,r=a.length;t<r;t++){D=a[t];if("M"==D[0])n=+D[1],p=+D[2];else{q=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6]);if(c+q>h){if(d&&!s.start){n=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6],h-c);O+=["C"+e(n.start.x),e(n.start.y),e(n.m.x),e(n.m.y),e(n.x),e(n.y)];if(l)return O;s.start=O;O=["M"+e(n.x),e(n.y)+"C"+e(n.n.x),e(n.n.y),e(n.end.x),e(n.end.y),e(D[5]),e(D[6])].join();c+=q;n=+D[5];p=+D[6];continue}if(!b&&!d)return n=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6],h-c)}c+=q;n=+D[5];p=+D[6]}O+=
D.shift()+D}s.end=O;return n=b?c:d?s:u(n,p,D[0],D[1],D[2],D[3],D[4],D[5],1)},null,a._.clone)}function u(a,b,d,e,h,f,k,l,n){var p=1-n,q=ma(p,3),s=ma(p,2),c=n*n,t=c*n,r=q*a+3*s*n*d+3*p*n*n*h+t*k,q=q*b+3*s*n*e+3*p*n*n*f+t*l,s=a+2*n*(d-a)+c*(h-2*d+a),t=b+2*n*(e-b)+c*(f-2*e+b),x=d+2*n*(h-d)+c*(k-2*h+d),c=e+2*n*(f-e)+c*(l-2*f+e);a=p*a+n*d;b=p*b+n*e;h=p*h+n*k;f=p*f+n*l;l=90-180*F.atan2(s-x,t-c)/S;return{x:r,y:q,m:{x:s,y:t},n:{x:x,y:c},start:{x:a,y:b},end:{x:h,y:f},alpha:l}}function p(b,d,e,h,f,n,k,l){a.is(b,
"array")||(b=[b,d,e,h,f,n,k,l]);b=U.apply(null,b);return w(b.min.x,b.min.y,b.max.x-b.min.x,b.max.y-b.min.y)}function b(a,b,d){return b>=a.x&&b<=a.x+a.width&&d>=a.y&&d<=a.y+a.height}function q(a,d){a=w(a);d=w(d);return b(d,a.x,a.y)||b(d,a.x2,a.y)||b(d,a.x,a.y2)||b(d,a.x2,a.y2)||b(a,d.x,d.y)||b(a,d.x2,d.y)||b(a,d.x,d.y2)||b(a,d.x2,d.y2)||(a.x<d.x2&&a.x>d.x||d.x<a.x2&&d.x>a.x)&&(a.y<d.y2&&a.y>d.y||d.y<a.y2&&d.y>a.y)}function e(a,b,d,e,h,f,n,k,l){null==l&&(l=1);l=(1<l?1:0>l?0:l)/2;for(var p=[-0.1252,
0.1252,-0.3678,0.3678,-0.5873,0.5873,-0.7699,0.7699,-0.9041,0.9041,-0.9816,0.9816],q=[0.2491,0.2491,0.2335,0.2335,0.2032,0.2032,0.1601,0.1601,0.1069,0.1069,0.0472,0.0472],s=0,c=0;12>c;c++)var t=l*p[c]+l,r=t*(t*(-3*a+9*d-9*h+3*n)+6*a-12*d+6*h)-3*a+3*d,t=t*(t*(-3*b+9*e-9*f+3*k)+6*b-12*e+6*f)-3*b+3*e,s=s+q[c]*F.sqrt(r*r+t*t);return l*s}function l(a,b,d){a=I(a);b=I(b);for(var h,f,l,n,k,s,r,O,x,c,t=d?0:[],w=0,v=a.length;w<v;w++)if(x=a[w],"M"==x[0])h=k=x[1],f=s=x[2];else{"C"==x[0]?(x=[h,f].concat(x.slice(1)),
h=x[6],f=x[7]):(x=[h,f,h,f,k,s,k,s],h=k,f=s);for(var G=0,y=b.length;G<y;G++)if(c=b[G],"M"==c[0])l=r=c[1],n=O=c[2];else{"C"==c[0]?(c=[l,n].concat(c.slice(1)),l=c[6],n=c[7]):(c=[l,n,l,n,r,O,r,O],l=r,n=O);var z;var K=x,B=c;z=d;var H=p(K),J=p(B);if(q(H,J)){for(var H=e.apply(0,K),J=e.apply(0,B),H=~~(H/8),J=~~(J/8),U=[],A=[],F={},M=z?0:[],P=0;P<H+1;P++){var C=u.apply(0,K.concat(P/H));U.push({x:C.x,y:C.y,t:P/H})}for(P=0;P<J+1;P++)C=u.apply(0,B.concat(P/J)),A.push({x:C.x,y:C.y,t:P/J});for(P=0;P<H;P++)for(K=
0;K<J;K++){var Q=U[P],L=U[P+1],B=A[K],C=A[K+1],N=0.001>Z(L.x-Q.x)?"y":"x",S=0.001>Z(C.x-B.x)?"y":"x",R;R=Q.x;var Y=Q.y,V=L.x,ea=L.y,fa=B.x,ga=B.y,ha=C.x,ia=C.y;if(W(R,V)<X(fa,ha)||X(R,V)>W(fa,ha)||W(Y,ea)<X(ga,ia)||X(Y,ea)>W(ga,ia))R=void 0;else{var $=(R*ea-Y*V)*(fa-ha)-(R-V)*(fa*ia-ga*ha),aa=(R*ea-Y*V)*(ga-ia)-(Y-ea)*(fa*ia-ga*ha),ja=(R-V)*(ga-ia)-(Y-ea)*(fa-ha);if(ja){var $=$/ja,aa=aa/ja,ja=+$.toFixed(2),ba=+aa.toFixed(2);R=ja<+X(R,V).toFixed(2)||ja>+W(R,V).toFixed(2)||ja<+X(fa,ha).toFixed(2)||
ja>+W(fa,ha).toFixed(2)||ba<+X(Y,ea).toFixed(2)||ba>+W(Y,ea).toFixed(2)||ba<+X(ga,ia).toFixed(2)||ba>+W(ga,ia).toFixed(2)?void 0:{x:$,y:aa}}else R=void 0}R&&F[R.x.toFixed(4)]!=R.y.toFixed(4)&&(F[R.x.toFixed(4)]=R.y.toFixed(4),Q=Q.t+Z((R[N]-Q[N])/(L[N]-Q[N]))*(L.t-Q.t),B=B.t+Z((R[S]-B[S])/(C[S]-B[S]))*(C.t-B.t),0<=Q&&1>=Q&&0<=B&&1>=B&&(z?M++:M.push({x:R.x,y:R.y,t1:Q,t2:B})))}z=M}else z=z?0:[];if(d)t+=z;else{H=0;for(J=z.length;H<J;H++)z[H].segment1=w,z[H].segment2=G,z[H].bez1=x,z[H].bez2=c;t=t.concat(z)}}}return t}
function r(a){var b=A(a);if(b.bbox)return C(b.bbox);if(!a)return w();a=I(a);for(var d=0,e=0,h=[],f=[],l,n=0,k=a.length;n<k;n++)l=a[n],"M"==l[0]?(d=l[1],e=l[2],h.push(d),f.push(e)):(d=U(d,e,l[1],l[2],l[3],l[4],l[5],l[6]),h=h.concat(d.min.x,d.max.x),f=f.concat(d.min.y,d.max.y),d=l[5],e=l[6]);a=X.apply(0,h);l=X.apply(0,f);h=W.apply(0,h);f=W.apply(0,f);f=w(a,l,h-a,f-l);b.bbox=C(f);return f}function s(a,b,d,e,h){if(h)return[["M",+a+ +h,b],["l",d-2*h,0],["a",h,h,0,0,1,h,h],["l",0,e-2*h],["a",h,h,0,0,1,
-h,h],["l",2*h-d,0],["a",h,h,0,0,1,-h,-h],["l",0,2*h-e],["a",h,h,0,0,1,h,-h],["z"] ];a=[["M",a,b],["l",d,0],["l",0,e],["l",-d,0],["z"] ];a.toString=z;return a}function x(a,b,d,e,h){null==h&&null==e&&(e=d);a=+a;b=+b;d=+d;e=+e;if(null!=h){var f=Math.PI/180,l=a+d*Math.cos(-e*f);a+=d*Math.cos(-h*f);var n=b+d*Math.sin(-e*f);b+=d*Math.sin(-h*f);d=[["M",l,n],["A",d,d,0,+(180<h-e),0,a,b] ]}else d=[["M",a,b],["m",0,-e],["a",d,e,0,1,1,0,2*e],["a",d,e,0,1,1,0,-2*e],["z"] ];d.toString=z;return d}function G(b){var e=
A(b);if(e.abs)return d(e.abs);Q(b,"array")&&Q(b&&b[0],"array")||(b=a.parsePathString(b));if(!b||!b.length)return[["M",0,0] ];var h=[],f=0,l=0,n=0,k=0,p=0;"M"==b[0][0]&&(f=+b[0][1],l=+b[0][2],n=f,k=l,p++,h[0]=["M",f,l]);for(var q=3==b.length&&"M"==b[0][0]&&"R"==b[1][0].toUpperCase()&&"Z"==b[2][0].toUpperCase(),s,r,w=p,c=b.length;w<c;w++){h.push(s=[]);r=b[w];p=r[0];if(p!=p.toUpperCase())switch(s[0]=p.toUpperCase(),s[0]){case "A":s[1]=r[1];s[2]=r[2];s[3]=r[3];s[4]=r[4];s[5]=r[5];s[6]=+r[6]+f;s[7]=+r[7]+
l;break;case "V":s[1]=+r[1]+l;break;case "H":s[1]=+r[1]+f;break;case "R":for(var t=[f,l].concat(r.slice(1)),u=2,v=t.length;u<v;u++)t[u]=+t[u]+f,t[++u]=+t[u]+l;h.pop();h=h.concat(P(t,q));break;case "O":h.pop();t=x(f,l,r[1],r[2]);t.push(t[0]);h=h.concat(t);break;case "U":h.pop();h=h.concat(x(f,l,r[1],r[2],r[3]));s=["U"].concat(h[h.length-1].slice(-2));break;case "M":n=+r[1]+f,k=+r[2]+l;default:for(u=1,v=r.length;u<v;u++)s[u]=+r[u]+(u%2?f:l)}else if("R"==p)t=[f,l].concat(r.slice(1)),h.pop(),h=h.concat(P(t,
q)),s=["R"].concat(r.slice(-2));else if("O"==p)h.pop(),t=x(f,l,r[1],r[2]),t.push(t[0]),h=h.concat(t);else if("U"==p)h.pop(),h=h.concat(x(f,l,r[1],r[2],r[3])),s=["U"].concat(h[h.length-1].slice(-2));else for(t=0,u=r.length;t<u;t++)s[t]=r[t];p=p.toUpperCase();if("O"!=p)switch(s[0]){case "Z":f=+n;l=+k;break;case "H":f=s[1];break;case "V":l=s[1];break;case "M":n=s[s.length-2],k=s[s.length-1];default:f=s[s.length-2],l=s[s.length-1]}}h.toString=z;e.abs=d(h);return h}function h(a,b,d,e){return[a,b,d,e,d,
e]}function J(a,b,d,e,h,f){var l=1/3,n=2/3;return[l*a+n*d,l*b+n*e,l*h+n*d,l*f+n*e,h,f]}function K(b,d,e,h,f,l,n,k,p,s){var r=120*S/180,q=S/180*(+f||0),c=[],t,x=a._.cacher(function(a,b,c){var d=a*F.cos(c)-b*F.sin(c);a=a*F.sin(c)+b*F.cos(c);return{x:d,y:a}});if(s)v=s[0],t=s[1],l=s[2],u=s[3];else{t=x(b,d,-q);b=t.x;d=t.y;t=x(k,p,-q);k=t.x;p=t.y;F.cos(S/180*f);F.sin(S/180*f);t=(b-k)/2;v=(d-p)/2;u=t*t/(e*e)+v*v/(h*h);1<u&&(u=F.sqrt(u),e*=u,h*=u);var u=e*e,w=h*h,u=(l==n?-1:1)*F.sqrt(Z((u*w-u*v*v-w*t*t)/
(u*v*v+w*t*t)));l=u*e*v/h+(b+k)/2;var u=u*-h*t/e+(d+p)/2,v=F.asin(((d-u)/h).toFixed(9));t=F.asin(((p-u)/h).toFixed(9));v=b<l?S-v:v;t=k<l?S-t:t;0>v&&(v=2*S+v);0>t&&(t=2*S+t);n&&v>t&&(v-=2*S);!n&&t>v&&(t-=2*S)}if(Z(t-v)>r){var c=t,w=k,G=p;t=v+r*(n&&t>v?1:-1);k=l+e*F.cos(t);p=u+h*F.sin(t);c=K(k,p,e,h,f,0,n,w,G,[t,c,l,u])}l=t-v;f=F.cos(v);r=F.sin(v);n=F.cos(t);t=F.sin(t);l=F.tan(l/4);e=4/3*e*l;l*=4/3*h;h=[b,d];b=[b+e*r,d-l*f];d=[k+e*t,p-l*n];k=[k,p];b[0]=2*h[0]-b[0];b[1]=2*h[1]-b[1];if(s)return[b,d,k].concat(c);
c=[b,d,k].concat(c).join().split(",");s=[];k=0;for(p=c.length;k<p;k++)s[k]=k%2?x(c[k-1],c[k],q).y:x(c[k],c[k+1],q).x;return s}function U(a,b,d,e,h,f,l,k){for(var n=[],p=[[],[] ],s,r,c,t,q=0;2>q;++q)0==q?(r=6*a-12*d+6*h,s=-3*a+9*d-9*h+3*l,c=3*d-3*a):(r=6*b-12*e+6*f,s=-3*b+9*e-9*f+3*k,c=3*e-3*b),1E-12>Z(s)?1E-12>Z(r)||(s=-c/r,0<s&&1>s&&n.push(s)):(t=r*r-4*c*s,c=F.sqrt(t),0>t||(t=(-r+c)/(2*s),0<t&&1>t&&n.push(t),s=(-r-c)/(2*s),0<s&&1>s&&n.push(s)));for(r=q=n.length;q--;)s=n[q],c=1-s,p[0][q]=c*c*c*a+3*
c*c*s*d+3*c*s*s*h+s*s*s*l,p[1][q]=c*c*c*b+3*c*c*s*e+3*c*s*s*f+s*s*s*k;p[0][r]=a;p[1][r]=b;p[0][r+1]=l;p[1][r+1]=k;p[0].length=p[1].length=r+2;return{min:{x:X.apply(0,p[0]),y:X.apply(0,p[1])},max:{x:W.apply(0,p[0]),y:W.apply(0,p[1])}}}function I(a,b){var e=!b&&A(a);if(!b&&e.curve)return d(e.curve);var f=G(a),l=b&&G(b),n={x:0,y:0,bx:0,by:0,X:0,Y:0,qx:null,qy:null},k={x:0,y:0,bx:0,by:0,X:0,Y:0,qx:null,qy:null},p=function(a,b,c){if(!a)return["C",b.x,b.y,b.x,b.y,b.x,b.y];a[0]in{T:1,Q:1}||(b.qx=b.qy=null);
switch(a[0]){case "M":b.X=a[1];b.Y=a[2];break;case "A":a=["C"].concat(K.apply(0,[b.x,b.y].concat(a.slice(1))));break;case "S":"C"==c||"S"==c?(c=2*b.x-b.bx,b=2*b.y-b.by):(c=b.x,b=b.y);a=["C",c,b].concat(a.slice(1));break;case "T":"Q"==c||"T"==c?(b.qx=2*b.x-b.qx,b.qy=2*b.y-b.qy):(b.qx=b.x,b.qy=b.y);a=["C"].concat(J(b.x,b.y,b.qx,b.qy,a[1],a[2]));break;case "Q":b.qx=a[1];b.qy=a[2];a=["C"].concat(J(b.x,b.y,a[1],a[2],a[3],a[4]));break;case "L":a=["C"].concat(h(b.x,b.y,a[1],a[2]));break;case "H":a=["C"].concat(h(b.x,
b.y,a[1],b.y));break;case "V":a=["C"].concat(h(b.x,b.y,b.x,a[1]));break;case "Z":a=["C"].concat(h(b.x,b.y,b.X,b.Y))}return a},s=function(a,b){if(7<a[b].length){a[b].shift();for(var c=a[b];c.length;)q[b]="A",l&&(u[b]="A"),a.splice(b++,0,["C"].concat(c.splice(0,6)));a.splice(b,1);v=W(f.length,l&&l.length||0)}},r=function(a,b,c,d,e){a&&b&&"M"==a[e][0]&&"M"!=b[e][0]&&(b.splice(e,0,["M",d.x,d.y]),c.bx=0,c.by=0,c.x=a[e][1],c.y=a[e][2],v=W(f.length,l&&l.length||0))},q=[],u=[],c="",t="",x=0,v=W(f.length,
l&&l.length||0);for(;x<v;x++){f[x]&&(c=f[x][0]);"C"!=c&&(q[x]=c,x&&(t=q[x-1]));f[x]=p(f[x],n,t);"A"!=q[x]&&"C"==c&&(q[x]="C");s(f,x);l&&(l[x]&&(c=l[x][0]),"C"!=c&&(u[x]=c,x&&(t=u[x-1])),l[x]=p(l[x],k,t),"A"!=u[x]&&"C"==c&&(u[x]="C"),s(l,x));r(f,l,n,k,x);r(l,f,k,n,x);var w=f[x],z=l&&l[x],y=w.length,U=l&&z.length;n.x=w[y-2];n.y=w[y-1];n.bx=$(w[y-4])||n.x;n.by=$(w[y-3])||n.y;k.bx=l&&($(z[U-4])||k.x);k.by=l&&($(z[U-3])||k.y);k.x=l&&z[U-2];k.y=l&&z[U-1]}l||(e.curve=d(f));return l?[f,l]:f}function P(a,
b){for(var d=[],e=0,h=a.length;h-2*!b>e;e+=2){var f=[{x:+a[e-2],y:+a[e-1]},{x:+a[e],y:+a[e+1]},{x:+a[e+2],y:+a[e+3]},{x:+a[e+4],y:+a[e+5]}];b?e?h-4==e?f[3]={x:+a[0],y:+a[1]}:h-2==e&&(f[2]={x:+a[0],y:+a[1]},f[3]={x:+a[2],y:+a[3]}):f[0]={x:+a[h-2],y:+a[h-1]}:h-4==e?f[3]=f[2]:e||(f[0]={x:+a[e],y:+a[e+1]});d.push(["C",(-f[0].x+6*f[1].x+f[2].x)/6,(-f[0].y+6*f[1].y+f[2].y)/6,(f[1].x+6*f[2].x-f[3].x)/6,(f[1].y+6*f[2].y-f[3].y)/6,f[2].x,f[2].y])}return d}y=k.prototype;var Q=a.is,C=a._.clone,L="hasOwnProperty",
N=/,?([a-z]),?/gi,$=parseFloat,F=Math,S=F.PI,X=F.min,W=F.max,ma=F.pow,Z=F.abs;M=n(1);var na=n(),ba=n(0,1),V=a._unit2px;a.path=A;a.path.getTotalLength=M;a.path.getPointAtLength=na;a.path.getSubpath=function(a,b,d){if(1E-6>this.getTotalLength(a)-d)return ba(a,b).end;a=ba(a,d,1);return b?ba(a,b).end:a};y.getTotalLength=function(){if(this.node.getTotalLength)return this.node.getTotalLength()};y.getPointAtLength=function(a){return na(this.attr("d"),a)};y.getSubpath=function(b,d){return a.path.getSubpath(this.attr("d"),
b,d)};a._.box=w;a.path.findDotsAtSegment=u;a.path.bezierBBox=p;a.path.isPointInsideBBox=b;a.path.isBBoxIntersect=q;a.path.intersection=function(a,b){return l(a,b)};a.path.intersectionNumber=function(a,b){return l(a,b,1)};a.path.isPointInside=function(a,d,e){var h=r(a);return b(h,d,e)&&1==l(a,[["M",d,e],["H",h.x2+10] ],1)%2};a.path.getBBox=r;a.path.get={path:function(a){return a.attr("path")},circle:function(a){a=V(a);return x(a.cx,a.cy,a.r)},ellipse:function(a){a=V(a);return x(a.cx||0,a.cy||0,a.rx,
a.ry)},rect:function(a){a=V(a);return s(a.x||0,a.y||0,a.width,a.height,a.rx,a.ry)},image:function(a){a=V(a);return s(a.x||0,a.y||0,a.width,a.height)},line:function(a){return"M"+[a.attr("x1")||0,a.attr("y1")||0,a.attr("x2"),a.attr("y2")]},polyline:function(a){return"M"+a.attr("points")},polygon:function(a){return"M"+a.attr("points")+"z"},deflt:function(a){a=a.node.getBBox();return s(a.x,a.y,a.width,a.height)}};a.path.toRelative=function(b){var e=A(b),h=String.prototype.toLowerCase;if(e.rel)return d(e.rel);
a.is(b,"array")&&a.is(b&&b[0],"array")||(b=a.parsePathString(b));var f=[],l=0,n=0,k=0,p=0,s=0;"M"==b[0][0]&&(l=b[0][1],n=b[0][2],k=l,p=n,s++,f.push(["M",l,n]));for(var r=b.length;s<r;s++){var q=f[s]=[],x=b[s];if(x[0]!=h.call(x[0]))switch(q[0]=h.call(x[0]),q[0]){case "a":q[1]=x[1];q[2]=x[2];q[3]=x[3];q[4]=x[4];q[5]=x[5];q[6]=+(x[6]-l).toFixed(3);q[7]=+(x[7]-n).toFixed(3);break;case "v":q[1]=+(x[1]-n).toFixed(3);break;case "m":k=x[1],p=x[2];default:for(var c=1,t=x.length;c<t;c++)q[c]=+(x[c]-(c%2?l:
n)).toFixed(3)}else for(f[s]=[],"m"==x[0]&&(k=x[1]+l,p=x[2]+n),q=0,c=x.length;q<c;q++)f[s][q]=x[q];x=f[s].length;switch(f[s][0]){case "z":l=k;n=p;break;case "h":l+=+f[s][x-1];break;case "v":n+=+f[s][x-1];break;default:l+=+f[s][x-2],n+=+f[s][x-1]}}f.toString=z;e.rel=d(f);return f};a.path.toAbsolute=G;a.path.toCubic=I;a.path.map=function(a,b){if(!b)return a;var d,e,h,f,l,n,k;a=I(a);h=0;for(l=a.length;h<l;h++)for(k=a[h],f=1,n=k.length;f<n;f+=2)d=b.x(k[f],k[f+1]),e=b.y(k[f],k[f+1]),k[f]=d,k[f+1]=e;return a};
a.path.toString=z;a.path.clone=d});C.plugin(function(a,v,y,C){var A=Math.max,w=Math.min,z=function(a){this.items=[];this.bindings={};this.length=0;this.type="set";if(a)for(var f=0,n=a.length;f<n;f++)a[f]&&(this[this.items.length]=this.items[this.items.length]=a[f],this.length++)};v=z.prototype;v.push=function(){for(var a,f,n=0,k=arguments.length;n<k;n++)if(a=arguments[n])f=this.items.length,this[f]=this.items[f]=a,this.length++;return this};v.pop=function(){this.length&&delete this[this.length--];
return this.items.pop()};v.forEach=function(a,f){for(var n=0,k=this.items.length;n<k&&!1!==a.call(f,this.items[n],n);n++);return this};v.animate=function(d,f,n,u){"function"!=typeof n||n.length||(u=n,n=L.linear);d instanceof a._.Animation&&(u=d.callback,n=d.easing,f=n.dur,d=d.attr);var p=arguments;if(a.is(d,"array")&&a.is(p[p.length-1],"array"))var b=!0;var q,e=function(){q?this.b=q:q=this.b},l=0,r=u&&function(){l++==this.length&&u.call(this)};return this.forEach(function(a,l){k.once("snap.animcreated."+
a.id,e);b?p[l]&&a.animate.apply(a,p[l]):a.animate(d,f,n,r)})};v.remove=function(){for(;this.length;)this.pop().remove();return this};v.bind=function(a,f,k){var u={};if("function"==typeof f)this.bindings[a]=f;else{var p=k||a;this.bindings[a]=function(a){u[p]=a;f.attr(u)}}return this};v.attr=function(a){var f={},k;for(k in a)if(this.bindings[k])this.bindings[k](a[k]);else f[k]=a[k];a=0;for(k=this.items.length;a<k;a++)this.items[a].attr(f);return this};v.clear=function(){for(;this.length;)this.pop()};
v.splice=function(a,f,k){a=0>a?A(this.length+a,0):a;f=A(0,w(this.length-a,f));var u=[],p=[],b=[],q;for(q=2;q<arguments.length;q++)b.push(arguments[q]);for(q=0;q<f;q++)p.push(this[a+q]);for(;q<this.length-a;q++)u.push(this[a+q]);var e=b.length;for(q=0;q<e+u.length;q++)this.items[a+q]=this[a+q]=q<e?b[q]:u[q-e];for(q=this.items.length=this.length-=f-e;this[q];)delete this[q++];return new z(p)};v.exclude=function(a){for(var f=0,k=this.length;f<k;f++)if(this[f]==a)return this.splice(f,1),!0;return!1};
v.insertAfter=function(a){for(var f=this.items.length;f--;)this.items[f].insertAfter(a);return this};v.getBBox=function(){for(var a=[],f=[],k=[],u=[],p=this.items.length;p--;)if(!this.items[p].removed){var b=this.items[p].getBBox();a.push(b.x);f.push(b.y);k.push(b.x+b.width);u.push(b.y+b.height)}a=w.apply(0,a);f=w.apply(0,f);k=A.apply(0,k);u=A.apply(0,u);return{x:a,y:f,x2:k,y2:u,width:k-a,height:u-f,cx:a+(k-a)/2,cy:f+(u-f)/2}};v.clone=function(a){a=new z;for(var f=0,k=this.items.length;f<k;f++)a.push(this.items[f].clone());
return a};v.toString=function(){return"Snap\u2018s set"};v.type="set";a.set=function(){var a=new z;arguments.length&&a.push.apply(a,Array.prototype.slice.call(arguments,0));return a}});C.plugin(function(a,v,y,C){function A(a){var b=a[0];switch(b.toLowerCase()){case "t":return[b,0,0];case "m":return[b,1,0,0,1,0,0];case "r":return 4==a.length?[b,0,a[2],a[3] ]:[b,0];case "s":return 5==a.length?[b,1,1,a[3],a[4] ]:3==a.length?[b,1,1]:[b,1]}}function w(b,d,f){d=q(d).replace(/\.{3}|\u2026/g,b);b=a.parseTransformString(b)||
[];d=a.parseTransformString(d)||[];for(var k=Math.max(b.length,d.length),p=[],v=[],h=0,w,z,y,I;h<k;h++){y=b[h]||A(d[h]);I=d[h]||A(y);if(y[0]!=I[0]||"r"==y[0].toLowerCase()&&(y[2]!=I[2]||y[3]!=I[3])||"s"==y[0].toLowerCase()&&(y[3]!=I[3]||y[4]!=I[4])){b=a._.transform2matrix(b,f());d=a._.transform2matrix(d,f());p=[["m",b.a,b.b,b.c,b.d,b.e,b.f] ];v=[["m",d.a,d.b,d.c,d.d,d.e,d.f] ];break}p[h]=[];v[h]=[];w=0;for(z=Math.max(y.length,I.length);w<z;w++)w in y&&(p[h][w]=y[w]),w in I&&(v[h][w]=I[w])}return{from:u(p),
to:u(v),f:n(p)}}function z(a){return a}function d(a){return function(b){return+b.toFixed(3)+a}}function f(b){return a.rgb(b[0],b[1],b[2])}function n(a){var b=0,d,f,k,n,h,p,q=[];d=0;for(f=a.length;d<f;d++){h="[";p=['"'+a[d][0]+'"'];k=1;for(n=a[d].length;k<n;k++)p[k]="val["+b++ +"]";h+=p+"]";q[d]=h}return Function("val","return Snap.path.toString.call(["+q+"])")}function u(a){for(var b=[],d=0,f=a.length;d<f;d++)for(var k=1,n=a[d].length;k<n;k++)b.push(a[d][k]);return b}var p={},b=/[a-z]+$/i,q=String;
p.stroke=p.fill="colour";v.prototype.equal=function(a,b){return k("snap.util.equal",this,a,b).firstDefined()};k.on("snap.util.equal",function(e,k){var r,s;r=q(this.attr(e)||"");var x=this;if(r==+r&&k==+k)return{from:+r,to:+k,f:z};if("colour"==p[e])return r=a.color(r),s=a.color(k),{from:[r.r,r.g,r.b,r.opacity],to:[s.r,s.g,s.b,s.opacity],f:f};if("transform"==e||"gradientTransform"==e||"patternTransform"==e)return k instanceof a.Matrix&&(k=k.toTransformString()),a._.rgTransform.test(k)||(k=a._.svgTransform2string(k)),
w(r,k,function(){return x.getBBox(1)});if("d"==e||"path"==e)return r=a.path.toCubic(r,k),{from:u(r[0]),to:u(r[1]),f:n(r[0])};if("points"==e)return r=q(r).split(a._.separator),s=q(k).split(a._.separator),{from:r,to:s,f:function(a){return a}};aUnit=r.match(b);s=q(k).match(b);return aUnit&&aUnit==s?{from:parseFloat(r),to:parseFloat(k),f:d(aUnit)}:{from:this.asPX(e),to:this.asPX(e,k),f:z}})});C.plugin(function(a,v,y,C){var A=v.prototype,w="createTouch"in C.doc;v="click dblclick mousedown mousemove mouseout mouseover mouseup touchstart touchmove touchend touchcancel".split(" ");
var z={mousedown:"touchstart",mousemove:"touchmove",mouseup:"touchend"},d=function(a,b){var d="y"==a?"scrollTop":"scrollLeft",e=b&&b.node?b.node.ownerDocument:C.doc;return e[d in e.documentElement?"documentElement":"body"][d]},f=function(){this.returnValue=!1},n=function(){return this.originalEvent.preventDefault()},u=function(){this.cancelBubble=!0},p=function(){return this.originalEvent.stopPropagation()},b=function(){if(C.doc.addEventListener)return function(a,b,e,f){var k=w&&z[b]?z[b]:b,l=function(k){var l=
d("y",f),q=d("x",f);if(w&&z.hasOwnProperty(b))for(var r=0,u=k.targetTouches&&k.targetTouches.length;r<u;r++)if(k.targetTouches[r].target==a||a.contains(k.targetTouches[r].target)){u=k;k=k.targetTouches[r];k.originalEvent=u;k.preventDefault=n;k.stopPropagation=p;break}return e.call(f,k,k.clientX+q,k.clientY+l)};b!==k&&a.addEventListener(b,l,!1);a.addEventListener(k,l,!1);return function(){b!==k&&a.removeEventListener(b,l,!1);a.removeEventListener(k,l,!1);return!0}};if(C.doc.attachEvent)return function(a,
b,e,h){var k=function(a){a=a||h.node.ownerDocument.window.event;var b=d("y",h),k=d("x",h),k=a.clientX+k,b=a.clientY+b;a.preventDefault=a.preventDefault||f;a.stopPropagation=a.stopPropagation||u;return e.call(h,a,k,b)};a.attachEvent("on"+b,k);return function(){a.detachEvent("on"+b,k);return!0}}}(),q=[],e=function(a){for(var b=a.clientX,e=a.clientY,f=d("y"),l=d("x"),n,p=q.length;p--;){n=q[p];if(w)for(var r=a.touches&&a.touches.length,u;r--;){if(u=a.touches[r],u.identifier==n.el._drag.id||n.el.node.contains(u.target)){b=
u.clientX;e=u.clientY;(a.originalEvent?a.originalEvent:a).preventDefault();break}}else a.preventDefault();b+=l;e+=f;k("snap.drag.move."+n.el.id,n.move_scope||n.el,b-n.el._drag.x,e-n.el._drag.y,b,e,a)}},l=function(b){a.unmousemove(e).unmouseup(l);for(var d=q.length,f;d--;)f=q[d],f.el._drag={},k("snap.drag.end."+f.el.id,f.end_scope||f.start_scope||f.move_scope||f.el,b);q=[]};for(y=v.length;y--;)(function(d){a[d]=A[d]=function(e,f){a.is(e,"function")&&(this.events=this.events||[],this.events.push({name:d,
f:e,unbind:b(this.node||document,d,e,f||this)}));return this};a["un"+d]=A["un"+d]=function(a){for(var b=this.events||[],e=b.length;e--;)if(b[e].name==d&&(b[e].f==a||!a)){b[e].unbind();b.splice(e,1);!b.length&&delete this.events;break}return this}})(v[y]);A.hover=function(a,b,d,e){return this.mouseover(a,d).mouseout(b,e||d)};A.unhover=function(a,b){return this.unmouseover(a).unmouseout(b)};var r=[];A.drag=function(b,d,f,h,n,p){function u(r,v,w){(r.originalEvent||r).preventDefault();this._drag.x=v;
this._drag.y=w;this._drag.id=r.identifier;!q.length&&a.mousemove(e).mouseup(l);q.push({el:this,move_scope:h,start_scope:n,end_scope:p});d&&k.on("snap.drag.start."+this.id,d);b&&k.on("snap.drag.move."+this.id,b);f&&k.on("snap.drag.end."+this.id,f);k("snap.drag.start."+this.id,n||h||this,v,w,r)}if(!arguments.length){var v;return this.drag(function(a,b){this.attr({transform:v+(v?"T":"t")+[a,b]})},function(){v=this.transform().local})}this._drag={};r.push({el:this,start:u});this.mousedown(u);return this};
A.undrag=function(){for(var b=r.length;b--;)r[b].el==this&&(this.unmousedown(r[b].start),r.splice(b,1),k.unbind("snap.drag.*."+this.id));!r.length&&a.unmousemove(e).unmouseup(l);return this}});C.plugin(function(a,v,y,C){y=y.prototype;var A=/^\s*url\((.+)\)/,w=String,z=a._.$;a.filter={};y.filter=function(d){var f=this;"svg"!=f.type&&(f=f.paper);d=a.parse(w(d));var k=a._.id(),u=z("filter");z(u,{id:k,filterUnits:"userSpaceOnUse"});u.appendChild(d.node);f.defs.appendChild(u);return new v(u)};k.on("snap.util.getattr.filter",
function(){k.stop();var d=z(this.node,"filter");if(d)return(d=w(d).match(A))&&a.select(d[1])});k.on("snap.util.attr.filter",function(d){if(d instanceof v&&"filter"==d.type){k.stop();var f=d.node.id;f||(z(d.node,{id:d.id}),f=d.id);z(this.node,{filter:a.url(f)})}d&&"none"!=d||(k.stop(),this.node.removeAttribute("filter"))});a.filter.blur=function(d,f){null==d&&(d=2);return a.format('<feGaussianBlur stdDeviation="{def}"/>',{def:null==f?d:[d,f]})};a.filter.blur.toString=function(){return this()};a.filter.shadow=
function(d,f,k,u,p){"string"==typeof k&&(p=u=k,k=4);"string"!=typeof u&&(p=u,u="#000");null==k&&(k=4);null==p&&(p=1);null==d&&(d=0,f=2);null==f&&(f=d);u=a.color(u||"#000");return a.format('<feGaussianBlur in="SourceAlpha" stdDeviation="{blur}"/><feOffset dx="{dx}" dy="{dy}" result="offsetblur"/><feFlood flood-color="{color}"/><feComposite in2="offsetblur" operator="in"/><feComponentTransfer><feFuncA type="linear" slope="{opacity}"/></feComponentTransfer><feMerge><feMergeNode/><feMergeNode in="SourceGraphic"/></feMerge>',
{color:u,dx:d,dy:f,blur:k,opacity:p})};a.filter.shadow.toString=function(){return this()};a.filter.grayscale=function(d){null==d&&(d=1);return a.format('<feColorMatrix type="matrix" values="{a} {b} {c} 0 0 {d} {e} {f} 0 0 {g} {b} {h} 0 0 0 0 0 1 0"/>',{a:0.2126+0.7874*(1-d),b:0.7152-0.7152*(1-d),c:0.0722-0.0722*(1-d),d:0.2126-0.2126*(1-d),e:0.7152+0.2848*(1-d),f:0.0722-0.0722*(1-d),g:0.2126-0.2126*(1-d),h:0.0722+0.9278*(1-d)})};a.filter.grayscale.toString=function(){return this()};a.filter.sepia=
function(d){null==d&&(d=1);return a.format('<feColorMatrix type="matrix" values="{a} {b} {c} 0 0 {d} {e} {f} 0 0 {g} {h} {i} 0 0 0 0 0 1 0"/>',{a:0.393+0.607*(1-d),b:0.769-0.769*(1-d),c:0.189-0.189*(1-d),d:0.349-0.349*(1-d),e:0.686+0.314*(1-d),f:0.168-0.168*(1-d),g:0.272-0.272*(1-d),h:0.534-0.534*(1-d),i:0.131+0.869*(1-d)})};a.filter.sepia.toString=function(){return this()};a.filter.saturate=function(d){null==d&&(d=1);return a.format('<feColorMatrix type="saturate" values="{amount}"/>',{amount:1-
d})};a.filter.saturate.toString=function(){return this()};a.filter.hueRotate=function(d){return a.format('<feColorMatrix type="hueRotate" values="{angle}"/>',{angle:d||0})};a.filter.hueRotate.toString=function(){return this()};a.filter.invert=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="table" tableValues="{amount} {amount2}"/><feFuncG type="table" tableValues="{amount} {amount2}"/><feFuncB type="table" tableValues="{amount} {amount2}"/></feComponentTransfer>',{amount:d,
amount2:1-d})};a.filter.invert.toString=function(){return this()};a.filter.brightness=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="linear" slope="{amount}"/><feFuncG type="linear" slope="{amount}"/><feFuncB type="linear" slope="{amount}"/></feComponentTransfer>',{amount:d})};a.filter.brightness.toString=function(){return this()};a.filter.contrast=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="linear" slope="{amount}" intercept="{amount2}"/><feFuncG type="linear" slope="{amount}" intercept="{amount2}"/><feFuncB type="linear" slope="{amount}" intercept="{amount2}"/></feComponentTransfer>',
{amount:d,amount2:0.5-d/2})};a.filter.contrast.toString=function(){return this()}});return C});

]]> </script>
</svg>



<div class="markdown"><h3>The <code>mat</code> interface</h3>
</div>


<div class="markdown"><p>The only &quot;must have&quot; interface is <code>mat</code>.</p>
</div>

<pre class='language-julia'><code class='language-julia'>mat(H)</code></pre>
<pre id='var-hash432223' class='code-output documenter-example-output'>2×2 Matrix{ComplexF64}:
 0.707107+0.0im   0.707107+0.0im
 0.707107+0.0im  -0.707107+0.0im</pre>


<div class="markdown"><p><code>PutBlock</code> is a composite block, its matrix representation is:</p>
<p class="tex">$$\texttt&#123;put&#40;n, i&#61;&gt;G&#41;&#125; &#61; I^&#123;\otimes n-i&#125; \otimes G \otimes I^&#123;\otimes i-1&#125;$$</p>
</div>

<pre class='language-julia'><code class='language-julia'>vizcircuit(put(3, 1=&gt;H))</code></pre>
<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:xlink="http://www.w3.org/1999/xlink"
     version="1.2"
     width="42.5mm" height="56.25mm" viewBox="0 0 42.5 56.25"
     stroke="none"
     fill="#000000"
     stroke-width="0.3"
     font-size="3.88"

     id="img-aa899102">
<defs>
  <marker id="arrow" markerWidth="15" markerHeight="7" refX="5" refY="3.5" orient="auto" markerUnits="strokeWidth">
    <path d="M0,0 L15,3.5 L0,7 z" stroke="context-stroke" fill="context-stroke"/>
  </marker>
</defs>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-aa899102-1">
  <g transform="translate(23.12,6.25)">
    <path d="M-5,-5 L5,-5 5,5 -5,5  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-aa899102-2">
  <g transform="translate(23.12,43.75)">
    <path fill="none" d="M10.62,0 L-10.62,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-aa899102-3">
  <g transform="translate(23.12,25)">
    <path fill="none" d="M10.62,0 L-10.62,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-aa899102-4">
  <g transform="translate(33.75,6.25)">
    <path fill="none" d="M0,0 L0,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-aa899102-5">
  <g transform="translate(30.94,6.25)">
    <path fill="none" d="M-2.81,0 L2.81,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-aa899102-6">
  <g transform="translate(15.31,6.25)">
    <path fill="none" d="M-2.81,0 L2.81,0 " class="primitive"/>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="5.64" id="img-aa899102-7">
  <g transform="translate(23.12,6.25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">H</text>
    </g>
  </g>
</g>
<script> <![CDATA[
(function(N){var k=/[\.\/]/,L=/\s*,\s*/,C=function(a,d){return a-d},a,v,y={n:{}},M=function(){for(var a=0,d=this.length;a<d;a++)if("undefined"!=typeof this[a])return this[a]},A=function(){for(var a=this.length;--a;)if("undefined"!=typeof this[a])return this[a]},w=function(k,d){k=String(k);var f=v,n=Array.prototype.slice.call(arguments,2),u=w.listeners(k),p=0,b,q=[],e={},l=[],r=a;l.firstDefined=M;l.lastDefined=A;a=k;for(var s=v=0,x=u.length;s<x;s++)"zIndex"in u[s]&&(q.push(u[s].zIndex),0>u[s].zIndex&&
(e[u[s].zIndex]=u[s]));for(q.sort(C);0>q[p];)if(b=e[q[p++] ],l.push(b.apply(d,n)),v)return v=f,l;for(s=0;s<x;s++)if(b=u[s],"zIndex"in b)if(b.zIndex==q[p]){l.push(b.apply(d,n));if(v)break;do if(p++,(b=e[q[p] ])&&l.push(b.apply(d,n)),v)break;while(b)}else e[b.zIndex]=b;else if(l.push(b.apply(d,n)),v)break;v=f;a=r;return l};w._events=y;w.listeners=function(a){a=a.split(k);var d=y,f,n,u,p,b,q,e,l=[d],r=[];u=0;for(p=a.length;u<p;u++){e=[];b=0;for(q=l.length;b<q;b++)for(d=l[b].n,f=[d[a[u] ],d["*"] ],n=2;n--;)if(d=
f[n])e.push(d),r=r.concat(d.f||[]);l=e}return r};w.on=function(a,d){a=String(a);if("function"!=typeof d)return function(){};for(var f=a.split(L),n=0,u=f.length;n<u;n++)(function(a){a=a.split(k);for(var b=y,f,e=0,l=a.length;e<l;e++)b=b.n,b=b.hasOwnProperty(a[e])&&b[a[e] ]||(b[a[e] ]={n:{}});b.f=b.f||[];e=0;for(l=b.f.length;e<l;e++)if(b.f[e]==d){f=!0;break}!f&&b.f.push(d)})(f[n]);return function(a){+a==+a&&(d.zIndex=+a)}};w.f=function(a){var d=[].slice.call(arguments,1);return function(){w.apply(null,
[a,null].concat(d).concat([].slice.call(arguments,0)))}};w.stop=function(){v=1};w.nt=function(k){return k?(new RegExp("(?:\\.|\\/|^)"+k+"(?:\\.|\\/|$)")).test(a):a};w.nts=function(){return a.split(k)};w.off=w.unbind=function(a,d){if(a){var f=a.split(L);if(1<f.length)for(var n=0,u=f.length;n<u;n++)w.off(f[n],d);else{for(var f=a.split(k),p,b,q,e,l=[y],n=0,u=f.length;n<u;n++)for(e=0;e<l.length;e+=q.length-2){q=[e,1];p=l[e].n;if("*"!=f[n])p[f[n] ]&&q.push(p[f[n] ]);else for(b in p)p.hasOwnProperty(b)&&
q.push(p[b]);l.splice.apply(l,q)}n=0;for(u=l.length;n<u;n++)for(p=l[n];p.n;){if(d){if(p.f){e=0;for(f=p.f.length;e<f;e++)if(p.f[e]==d){p.f.splice(e,1);break}!p.f.length&&delete p.f}for(b in p.n)if(p.n.hasOwnProperty(b)&&p.n[b].f){q=p.n[b].f;e=0;for(f=q.length;e<f;e++)if(q[e]==d){q.splice(e,1);break}!q.length&&delete p.n[b].f}}else for(b in delete p.f,p.n)p.n.hasOwnProperty(b)&&p.n[b].f&&delete p.n[b].f;p=p.n}}}else w._events=y={n:{}}};w.once=function(a,d){var f=function(){w.unbind(a,f);return d.apply(this,
arguments)};return w.on(a,f)};w.version="0.4.2";w.toString=function(){return"You are running Eve 0.4.2"};"undefined"!=typeof module&&module.exports?module.exports=w:"function"===typeof define&&define.amd?define("eve",[],function(){return w}):N.eve=w})(this);
(function(N,k){"function"===typeof define&&define.amd?define("Snap.svg",["eve"],function(L){return k(N,L)}):k(N,N.eve)})(this,function(N,k){var L=function(a){var k={},y=N.requestAnimationFrame||N.webkitRequestAnimationFrame||N.mozRequestAnimationFrame||N.oRequestAnimationFrame||N.msRequestAnimationFrame||function(a){setTimeout(a,16)},M=Array.isArray||function(a){return a instanceof Array||"[object Array]"==Object.prototype.toString.call(a)},A=0,w="M"+(+new Date).toString(36),z=function(a){if(null==
a)return this.s;var b=this.s-a;this.b+=this.dur*b;this.B+=this.dur*b;this.s=a},d=function(a){if(null==a)return this.spd;this.spd=a},f=function(a){if(null==a)return this.dur;this.s=this.s*a/this.dur;this.dur=a},n=function(){delete k[this.id];this.update();a("mina.stop."+this.id,this)},u=function(){this.pdif||(delete k[this.id],this.update(),this.pdif=this.get()-this.b)},p=function(){this.pdif&&(this.b=this.get()-this.pdif,delete this.pdif,k[this.id]=this)},b=function(){var a;if(M(this.start)){a=[];
for(var b=0,e=this.start.length;b<e;b++)a[b]=+this.start[b]+(this.end[b]-this.start[b])*this.easing(this.s)}else a=+this.start+(this.end-this.start)*this.easing(this.s);this.set(a)},q=function(){var l=0,b;for(b in k)if(k.hasOwnProperty(b)){var e=k[b],f=e.get();l++;e.s=(f-e.b)/(e.dur/e.spd);1<=e.s&&(delete k[b],e.s=1,l--,function(b){setTimeout(function(){a("mina.finish."+b.id,b)})}(e));e.update()}l&&y(q)},e=function(a,r,s,x,G,h,J){a={id:w+(A++).toString(36),start:a,end:r,b:s,s:0,dur:x-s,spd:1,get:G,
set:h,easing:J||e.linear,status:z,speed:d,duration:f,stop:n,pause:u,resume:p,update:b};k[a.id]=a;r=0;for(var K in k)if(k.hasOwnProperty(K)&&(r++,2==r))break;1==r&&y(q);return a};e.time=Date.now||function(){return+new Date};e.getById=function(a){return k[a]||null};e.linear=function(a){return a};e.easeout=function(a){return Math.pow(a,1.7)};e.easein=function(a){return Math.pow(a,0.48)};e.easeinout=function(a){if(1==a)return 1;if(0==a)return 0;var b=0.48-a/1.04,e=Math.sqrt(0.1734+b*b);a=e-b;a=Math.pow(Math.abs(a),
1/3)*(0>a?-1:1);b=-e-b;b=Math.pow(Math.abs(b),1/3)*(0>b?-1:1);a=a+b+0.5;return 3*(1-a)*a*a+a*a*a};e.backin=function(a){return 1==a?1:a*a*(2.70158*a-1.70158)};e.backout=function(a){if(0==a)return 0;a-=1;return a*a*(2.70158*a+1.70158)+1};e.elastic=function(a){return a==!!a?a:Math.pow(2,-10*a)*Math.sin(2*(a-0.075)*Math.PI/0.3)+1};e.bounce=function(a){a<1/2.75?a*=7.5625*a:a<2/2.75?(a-=1.5/2.75,a=7.5625*a*a+0.75):a<2.5/2.75?(a-=2.25/2.75,a=7.5625*a*a+0.9375):(a-=2.625/2.75,a=7.5625*a*a+0.984375);return a};
return N.mina=e}("undefined"==typeof k?function(){}:k),C=function(){function a(c,t){if(c){if(c.tagName)return x(c);if(y(c,"array")&&a.set)return a.set.apply(a,c);if(c instanceof e)return c;if(null==t)return c=G.doc.querySelector(c),x(c)}return new s(null==c?"100%":c,null==t?"100%":t)}function v(c,a){if(a){"#text"==c&&(c=G.doc.createTextNode(a.text||""));"string"==typeof c&&(c=v(c));if("string"==typeof a)return"xlink:"==a.substring(0,6)?c.getAttributeNS(m,a.substring(6)):"xml:"==a.substring(0,4)?c.getAttributeNS(la,
a.substring(4)):c.getAttribute(a);for(var da in a)if(a[h](da)){var b=J(a[da]);b?"xlink:"==da.substring(0,6)?c.setAttributeNS(m,da.substring(6),b):"xml:"==da.substring(0,4)?c.setAttributeNS(la,da.substring(4),b):c.setAttribute(da,b):c.removeAttribute(da)}}else c=G.doc.createElementNS(la,c);return c}function y(c,a){a=J.prototype.toLowerCase.call(a);return"finite"==a?isFinite(c):"array"==a&&(c instanceof Array||Array.isArray&&Array.isArray(c))?!0:"null"==a&&null===c||a==typeof c&&null!==c||"object"==
a&&c===Object(c)||$.call(c).slice(8,-1).toLowerCase()==a}function M(c){if("function"==typeof c||Object(c)!==c)return c;var a=new c.constructor,b;for(b in c)c[h](b)&&(a[b]=M(c[b]));return a}function A(c,a,b){function m(){var e=Array.prototype.slice.call(arguments,0),f=e.join("\u2400"),d=m.cache=m.cache||{},l=m.count=m.count||[];if(d[h](f)){a:for(var e=l,l=f,B=0,H=e.length;B<H;B++)if(e[B]===l){e.push(e.splice(B,1)[0]);break a}return b?b(d[f]):d[f]}1E3<=l.length&&delete d[l.shift()];l.push(f);d[f]=c.apply(a,
e);return b?b(d[f]):d[f]}return m}function w(c,a,b,m,e,f){return null==e?(c-=b,a-=m,c||a?(180*I.atan2(-a,-c)/C+540)%360:0):w(c,a,e,f)-w(b,m,e,f)}function z(c){return c%360*C/180}function d(c){var a=[];c=c.replace(/(?:^|\s)(\w+)\(([^)]+)\)/g,function(c,b,m){m=m.split(/\s*,\s*|\s+/);"rotate"==b&&1==m.length&&m.push(0,0);"scale"==b&&(2<m.length?m=m.slice(0,2):2==m.length&&m.push(0,0),1==m.length&&m.push(m[0],0,0));"skewX"==b?a.push(["m",1,0,I.tan(z(m[0])),1,0,0]):"skewY"==b?a.push(["m",1,I.tan(z(m[0])),
0,1,0,0]):a.push([b.charAt(0)].concat(m));return c});return a}function f(c,t){var b=O(c),m=new a.Matrix;if(b)for(var e=0,f=b.length;e<f;e++){var h=b[e],d=h.length,B=J(h[0]).toLowerCase(),H=h[0]!=B,l=H?m.invert():0,E;"t"==B&&2==d?m.translate(h[1],0):"t"==B&&3==d?H?(d=l.x(0,0),B=l.y(0,0),H=l.x(h[1],h[2]),l=l.y(h[1],h[2]),m.translate(H-d,l-B)):m.translate(h[1],h[2]):"r"==B?2==d?(E=E||t,m.rotate(h[1],E.x+E.width/2,E.y+E.height/2)):4==d&&(H?(H=l.x(h[2],h[3]),l=l.y(h[2],h[3]),m.rotate(h[1],H,l)):m.rotate(h[1],
h[2],h[3])):"s"==B?2==d||3==d?(E=E||t,m.scale(h[1],h[d-1],E.x+E.width/2,E.y+E.height/2)):4==d?H?(H=l.x(h[2],h[3]),l=l.y(h[2],h[3]),m.scale(h[1],h[1],H,l)):m.scale(h[1],h[1],h[2],h[3]):5==d&&(H?(H=l.x(h[3],h[4]),l=l.y(h[3],h[4]),m.scale(h[1],h[2],H,l)):m.scale(h[1],h[2],h[3],h[4])):"m"==B&&7==d&&m.add(h[1],h[2],h[3],h[4],h[5],h[6])}return m}function n(c,t){if(null==t){var m=!0;t="linearGradient"==c.type||"radialGradient"==c.type?c.node.getAttribute("gradientTransform"):"pattern"==c.type?c.node.getAttribute("patternTransform"):
c.node.getAttribute("transform");if(!t)return new a.Matrix;t=d(t)}else t=a._.rgTransform.test(t)?J(t).replace(/\.{3}|\u2026/g,c._.transform||aa):d(t),y(t,"array")&&(t=a.path?a.path.toString.call(t):J(t)),c._.transform=t;var b=f(t,c.getBBox(1));if(m)return b;c.matrix=b}function u(c){c=c.node.ownerSVGElement&&x(c.node.ownerSVGElement)||c.node.parentNode&&x(c.node.parentNode)||a.select("svg")||a(0,0);var t=c.select("defs"),t=null==t?!1:t.node;t||(t=r("defs",c.node).node);return t}function p(c){return c.node.ownerSVGElement&&
x(c.node.ownerSVGElement)||a.select("svg")}function b(c,a,m){function b(c){if(null==c)return aa;if(c==+c)return c;v(B,{width:c});try{return B.getBBox().width}catch(a){return 0}}function h(c){if(null==c)return aa;if(c==+c)return c;v(B,{height:c});try{return B.getBBox().height}catch(a){return 0}}function e(b,B){null==a?d[b]=B(c.attr(b)||0):b==a&&(d=B(null==m?c.attr(b)||0:m))}var f=p(c).node,d={},B=f.querySelector(".svg---mgr");B||(B=v("rect"),v(B,{x:-9E9,y:-9E9,width:10,height:10,"class":"svg---mgr",
fill:"none"}),f.appendChild(B));switch(c.type){case "rect":e("rx",b),e("ry",h);case "image":e("width",b),e("height",h);case "text":e("x",b);e("y",h);break;case "circle":e("cx",b);e("cy",h);e("r",b);break;case "ellipse":e("cx",b);e("cy",h);e("rx",b);e("ry",h);break;case "line":e("x1",b);e("x2",b);e("y1",h);e("y2",h);break;case "marker":e("refX",b);e("markerWidth",b);e("refY",h);e("markerHeight",h);break;case "radialGradient":e("fx",b);e("fy",h);break;case "tspan":e("dx",b);e("dy",h);break;default:e(a,
b)}f.removeChild(B);return d}function q(c){y(c,"array")||(c=Array.prototype.slice.call(arguments,0));for(var a=0,b=0,m=this.node;this[a];)delete this[a++];for(a=0;a<c.length;a++)"set"==c[a].type?c[a].forEach(function(c){m.appendChild(c.node)}):m.appendChild(c[a].node);for(var h=m.childNodes,a=0;a<h.length;a++)this[b++]=x(h[a]);return this}function e(c){if(c.snap in E)return E[c.snap];var a=this.id=V(),b;try{b=c.ownerSVGElement}catch(m){}this.node=c;b&&(this.paper=new s(b));this.type=c.tagName;this.anims=
{};this._={transform:[]};c.snap=a;E[a]=this;"g"==this.type&&(this.add=q);if(this.type in{g:1,mask:1,pattern:1})for(var e in s.prototype)s.prototype[h](e)&&(this[e]=s.prototype[e])}function l(c){this.node=c}function r(c,a){var b=v(c);a.appendChild(b);return x(b)}function s(c,a){var b,m,f,d=s.prototype;if(c&&"svg"==c.tagName){if(c.snap in E)return E[c.snap];var l=c.ownerDocument;b=new e(c);m=c.getElementsByTagName("desc")[0];f=c.getElementsByTagName("defs")[0];m||(m=v("desc"),m.appendChild(l.createTextNode("Created with Snap")),
b.node.appendChild(m));f||(f=v("defs"),b.node.appendChild(f));b.defs=f;for(var ca in d)d[h](ca)&&(b[ca]=d[ca]);b.paper=b.root=b}else b=r("svg",G.doc.body),v(b.node,{height:a,version:1.1,width:c,xmlns:la});return b}function x(c){return!c||c instanceof e||c instanceof l?c:c.tagName&&"svg"==c.tagName.toLowerCase()?new s(c):c.tagName&&"object"==c.tagName.toLowerCase()&&"image/svg+xml"==c.type?new s(c.contentDocument.getElementsByTagName("svg")[0]):new e(c)}a.version="0.3.0";a.toString=function(){return"Snap v"+
this.version};a._={};var G={win:N,doc:N.document};a._.glob=G;var h="hasOwnProperty",J=String,K=parseFloat,U=parseInt,I=Math,P=I.max,Q=I.min,Y=I.abs,C=I.PI,aa="",$=Object.prototype.toString,F=/^\s*((#[a-f\d]{6})|(#[a-f\d]{3})|rgba?\(\s*([\d\.]+%?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+%?(?:\s*,\s*[\d\.]+%?)?)\s*\)|hsba?\(\s*([\d\.]+(?:deg|\xb0|%)?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+(?:%?\s*,\s*[\d\.]+)?%?)\s*\)|hsla?\(\s*([\d\.]+(?:deg|\xb0|%)?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+(?:%?\s*,\s*[\d\.]+)?%?)\s*\))\s*$/i;a._.separator=
RegExp("[,\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]+");var S=RegExp("[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*"),X={hs:1,rg:1},W=RegExp("([a-z])[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029,]*((-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*)+)",
"ig"),ma=RegExp("([rstm])[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029,]*((-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*)+)","ig"),Z=RegExp("(-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?)[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*",
"ig"),na=0,ba="S"+(+new Date).toString(36),V=function(){return ba+(na++).toString(36)},m="http://www.w3.org/1999/xlink",la="http://www.w3.org/2000/svg",E={},ca=a.url=function(c){return"url('#"+c+"')"};a._.$=v;a._.id=V;a.format=function(){var c=/\{([^\}]+)\}/g,a=/(?:(?:^|\.)(.+?)(?=\[|\.|$|\()|\[('|")(.+?)\2\])(\(\))?/g,b=function(c,b,m){var h=m;b.replace(a,function(c,a,b,m,t){a=a||m;h&&(a in h&&(h=h[a]),"function"==typeof h&&t&&(h=h()))});return h=(null==h||h==m?c:h)+""};return function(a,m){return J(a).replace(c,
function(c,a){return b(c,a,m)})}}();a._.clone=M;a._.cacher=A;a.rad=z;a.deg=function(c){return 180*c/C%360};a.angle=w;a.is=y;a.snapTo=function(c,a,b){b=y(b,"finite")?b:10;if(y(c,"array"))for(var m=c.length;m--;){if(Y(c[m]-a)<=b)return c[m]}else{c=+c;m=a%c;if(m<b)return a-m;if(m>c-b)return a-m+c}return a};a.getRGB=A(function(c){if(!c||(c=J(c)).indexOf("-")+1)return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka};if("none"==c)return{r:-1,g:-1,b:-1,hex:"none",toString:ka};!X[h](c.toLowerCase().substring(0,
2))&&"#"!=c.charAt()&&(c=T(c));if(!c)return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka};var b,m,e,f,d;if(c=c.match(F)){c[2]&&(e=U(c[2].substring(5),16),m=U(c[2].substring(3,5),16),b=U(c[2].substring(1,3),16));c[3]&&(e=U((d=c[3].charAt(3))+d,16),m=U((d=c[3].charAt(2))+d,16),b=U((d=c[3].charAt(1))+d,16));c[4]&&(d=c[4].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b*=2.55),m=K(d[1]),"%"==d[1].slice(-1)&&(m*=2.55),e=K(d[2]),"%"==d[2].slice(-1)&&(e*=2.55),"rgba"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),
d[3]&&"%"==d[3].slice(-1)&&(f/=100));if(c[5])return d=c[5].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b/=100),m=K(d[1]),"%"==d[1].slice(-1)&&(m/=100),e=K(d[2]),"%"==d[2].slice(-1)&&(e/=100),"deg"!=d[0].slice(-3)&&"\u00b0"!=d[0].slice(-1)||(b/=360),"hsba"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),d[3]&&"%"==d[3].slice(-1)&&(f/=100),a.hsb2rgb(b,m,e,f);if(c[6])return d=c[6].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b/=100),m=K(d[1]),"%"==d[1].slice(-1)&&(m/=100),e=K(d[2]),"%"==d[2].slice(-1)&&(e/=100),
"deg"!=d[0].slice(-3)&&"\u00b0"!=d[0].slice(-1)||(b/=360),"hsla"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),d[3]&&"%"==d[3].slice(-1)&&(f/=100),a.hsl2rgb(b,m,e,f);b=Q(I.round(b),255);m=Q(I.round(m),255);e=Q(I.round(e),255);f=Q(P(f,0),1);c={r:b,g:m,b:e,toString:ka};c.hex="#"+(16777216|e|m<<8|b<<16).toString(16).slice(1);c.opacity=y(f,"finite")?f:1;return c}return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka}},a);a.hsb=A(function(c,b,m){return a.hsb2rgb(c,b,m).hex});a.hsl=A(function(c,b,m){return a.hsl2rgb(c,
b,m).hex});a.rgb=A(function(c,a,b,m){if(y(m,"finite")){var e=I.round;return"rgba("+[e(c),e(a),e(b),+m.toFixed(2)]+")"}return"#"+(16777216|b|a<<8|c<<16).toString(16).slice(1)});var T=function(c){var a=G.doc.getElementsByTagName("head")[0]||G.doc.getElementsByTagName("svg")[0];T=A(function(c){if("red"==c.toLowerCase())return"rgb(255, 0, 0)";a.style.color="rgb(255, 0, 0)";a.style.color=c;c=G.doc.defaultView.getComputedStyle(a,aa).getPropertyValue("color");return"rgb(255, 0, 0)"==c?null:c});return T(c)},
qa=function(){return"hsb("+[this.h,this.s,this.b]+")"},ra=function(){return"hsl("+[this.h,this.s,this.l]+")"},ka=function(){return 1==this.opacity||null==this.opacity?this.hex:"rgba("+[this.r,this.g,this.b,this.opacity]+")"},D=function(c,b,m){null==b&&y(c,"object")&&"r"in c&&"g"in c&&"b"in c&&(m=c.b,b=c.g,c=c.r);null==b&&y(c,string)&&(m=a.getRGB(c),c=m.r,b=m.g,m=m.b);if(1<c||1<b||1<m)c/=255,b/=255,m/=255;return[c,b,m]},oa=function(c,b,m,e){c=I.round(255*c);b=I.round(255*b);m=I.round(255*m);c={r:c,
g:b,b:m,opacity:y(e,"finite")?e:1,hex:a.rgb(c,b,m),toString:ka};y(e,"finite")&&(c.opacity=e);return c};a.color=function(c){var b;y(c,"object")&&"h"in c&&"s"in c&&"b"in c?(b=a.hsb2rgb(c),c.r=b.r,c.g=b.g,c.b=b.b,c.opacity=1,c.hex=b.hex):y(c,"object")&&"h"in c&&"s"in c&&"l"in c?(b=a.hsl2rgb(c),c.r=b.r,c.g=b.g,c.b=b.b,c.opacity=1,c.hex=b.hex):(y(c,"string")&&(c=a.getRGB(c)),y(c,"object")&&"r"in c&&"g"in c&&"b"in c&&!("error"in c)?(b=a.rgb2hsl(c),c.h=b.h,c.s=b.s,c.l=b.l,b=a.rgb2hsb(c),c.v=b.b):(c={hex:"none"},
c.r=c.g=c.b=c.h=c.s=c.v=c.l=-1,c.error=1));c.toString=ka;return c};a.hsb2rgb=function(c,a,b,m){y(c,"object")&&"h"in c&&"s"in c&&"b"in c&&(b=c.b,a=c.s,c=c.h,m=c.o);var e,h,d;c=360*c%360/60;d=b*a;a=d*(1-Y(c%2-1));b=e=h=b-d;c=~~c;b+=[d,a,0,0,a,d][c];e+=[a,d,d,a,0,0][c];h+=[0,0,a,d,d,a][c];return oa(b,e,h,m)};a.hsl2rgb=function(c,a,b,m){y(c,"object")&&"h"in c&&"s"in c&&"l"in c&&(b=c.l,a=c.s,c=c.h);if(1<c||1<a||1<b)c/=360,a/=100,b/=100;var e,h,d;c=360*c%360/60;d=2*a*(0.5>b?b:1-b);a=d*(1-Y(c%2-1));b=e=
h=b-d/2;c=~~c;b+=[d,a,0,0,a,d][c];e+=[a,d,d,a,0,0][c];h+=[0,0,a,d,d,a][c];return oa(b,e,h,m)};a.rgb2hsb=function(c,a,b){b=D(c,a,b);c=b[0];a=b[1];b=b[2];var m,e;m=P(c,a,b);e=m-Q(c,a,b);c=((0==e?0:m==c?(a-b)/e:m==a?(b-c)/e+2:(c-a)/e+4)+360)%6*60/360;return{h:c,s:0==e?0:e/m,b:m,toString:qa}};a.rgb2hsl=function(c,a,b){b=D(c,a,b);c=b[0];a=b[1];b=b[2];var m,e,h;m=P(c,a,b);e=Q(c,a,b);h=m-e;c=((0==h?0:m==c?(a-b)/h:m==a?(b-c)/h+2:(c-a)/h+4)+360)%6*60/360;m=(m+e)/2;return{h:c,s:0==h?0:0.5>m?h/(2*m):h/(2-2*
m),l:m,toString:ra}};a.parsePathString=function(c){if(!c)return null;var b=a.path(c);if(b.arr)return a.path.clone(b.arr);var m={a:7,c:6,o:2,h:1,l:2,m:2,r:4,q:4,s:4,t:2,v:1,u:3,z:0},e=[];y(c,"array")&&y(c[0],"array")&&(e=a.path.clone(c));e.length||J(c).replace(W,function(c,a,b){var h=[];c=a.toLowerCase();b.replace(Z,function(c,a){a&&h.push(+a)});"m"==c&&2<h.length&&(e.push([a].concat(h.splice(0,2))),c="l",a="m"==a?"l":"L");"o"==c&&1==h.length&&e.push([a,h[0] ]);if("r"==c)e.push([a].concat(h));else for(;h.length>=
m[c]&&(e.push([a].concat(h.splice(0,m[c]))),m[c]););});e.toString=a.path.toString;b.arr=a.path.clone(e);return e};var O=a.parseTransformString=function(c){if(!c)return null;var b=[];y(c,"array")&&y(c[0],"array")&&(b=a.path.clone(c));b.length||J(c).replace(ma,function(c,a,m){var e=[];a.toLowerCase();m.replace(Z,function(c,a){a&&e.push(+a)});b.push([a].concat(e))});b.toString=a.path.toString;return b};a._.svgTransform2string=d;a._.rgTransform=RegExp("^[a-z][\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*-?\\.?\\d",
"i");a._.transform2matrix=f;a._unit2px=b;a._.getSomeDefs=u;a._.getSomeSVG=p;a.select=function(c){return x(G.doc.querySelector(c))};a.selectAll=function(c){c=G.doc.querySelectorAll(c);for(var b=(a.set||Array)(),m=0;m<c.length;m++)b.push(x(c[m]));return b};setInterval(function(){for(var c in E)if(E[h](c)){var a=E[c],b=a.node;("svg"!=a.type&&!b.ownerSVGElement||"svg"==a.type&&(!b.parentNode||"ownerSVGElement"in b.parentNode&&!b.ownerSVGElement))&&delete E[c]}},1E4);(function(c){function m(c){function a(c,
b){var m=v(c.node,b);(m=(m=m&&m.match(d))&&m[2])&&"#"==m.charAt()&&(m=m.substring(1))&&(f[m]=(f[m]||[]).concat(function(a){var m={};m[b]=ca(a);v(c.node,m)}))}function b(c){var a=v(c.node,"xlink:href");a&&"#"==a.charAt()&&(a=a.substring(1))&&(f[a]=(f[a]||[]).concat(function(a){c.attr("xlink:href","#"+a)}))}var e=c.selectAll("*"),h,d=/^\s*url\(("|'|)(.*)\1\)\s*$/;c=[];for(var f={},l=0,E=e.length;l<E;l++){h=e[l];a(h,"fill");a(h,"stroke");a(h,"filter");a(h,"mask");a(h,"clip-path");b(h);var t=v(h.node,
"id");t&&(v(h.node,{id:h.id}),c.push({old:t,id:h.id}))}l=0;for(E=c.length;l<E;l++)if(e=f[c[l].old])for(h=0,t=e.length;h<t;h++)e[h](c[l].id)}function e(c,a,b){return function(m){m=m.slice(c,a);1==m.length&&(m=m[0]);return b?b(m):m}}function d(c){return function(){var a=c?"<"+this.type:"",b=this.node.attributes,m=this.node.childNodes;if(c)for(var e=0,h=b.length;e<h;e++)a+=" "+b[e].name+'="'+b[e].value.replace(/"/g,'\\"')+'"';if(m.length){c&&(a+=">");e=0;for(h=m.length;e<h;e++)3==m[e].nodeType?a+=m[e].nodeValue:
1==m[e].nodeType&&(a+=x(m[e]).toString());c&&(a+="</"+this.type+">")}else c&&(a+="/>");return a}}c.attr=function(c,a){if(!c)return this;if(y(c,"string"))if(1<arguments.length){var b={};b[c]=a;c=b}else return k("snap.util.getattr."+c,this).firstDefined();for(var m in c)c[h](m)&&k("snap.util.attr."+m,this,c[m]);return this};c.getBBox=function(c){if(!a.Matrix||!a.path)return this.node.getBBox();var b=this,m=new a.Matrix;if(b.removed)return a._.box();for(;"use"==b.type;)if(c||(m=m.add(b.transform().localMatrix.translate(b.attr("x")||
0,b.attr("y")||0))),b.original)b=b.original;else var e=b.attr("xlink:href"),b=b.original=b.node.ownerDocument.getElementById(e.substring(e.indexOf("#")+1));var e=b._,h=a.path.get[b.type]||a.path.get.deflt;try{if(c)return e.bboxwt=h?a.path.getBBox(b.realPath=h(b)):a._.box(b.node.getBBox()),a._.box(e.bboxwt);b.realPath=h(b);b.matrix=b.transform().localMatrix;e.bbox=a.path.getBBox(a.path.map(b.realPath,m.add(b.matrix)));return a._.box(e.bbox)}catch(d){return a._.box()}};var f=function(){return this.string};
c.transform=function(c){var b=this._;if(null==c){var m=this;c=new a.Matrix(this.node.getCTM());for(var e=n(this),h=[e],d=new a.Matrix,l=e.toTransformString(),b=J(e)==J(this.matrix)?J(b.transform):l;"svg"!=m.type&&(m=m.parent());)h.push(n(m));for(m=h.length;m--;)d.add(h[m]);return{string:b,globalMatrix:c,totalMatrix:d,localMatrix:e,diffMatrix:c.clone().add(e.invert()),global:c.toTransformString(),total:d.toTransformString(),local:l,toString:f}}c instanceof a.Matrix?this.matrix=c:n(this,c);this.node&&
("linearGradient"==this.type||"radialGradient"==this.type?v(this.node,{gradientTransform:this.matrix}):"pattern"==this.type?v(this.node,{patternTransform:this.matrix}):v(this.node,{transform:this.matrix}));return this};c.parent=function(){return x(this.node.parentNode)};c.append=c.add=function(c){if(c){if("set"==c.type){var a=this;c.forEach(function(c){a.add(c)});return this}c=x(c);this.node.appendChild(c.node);c.paper=this.paper}return this};c.appendTo=function(c){c&&(c=x(c),c.append(this));return this};
c.prepend=function(c){if(c){if("set"==c.type){var a=this,b;c.forEach(function(c){b?b.after(c):a.prepend(c);b=c});return this}c=x(c);var m=c.parent();this.node.insertBefore(c.node,this.node.firstChild);this.add&&this.add();c.paper=this.paper;this.parent()&&this.parent().add();m&&m.add()}return this};c.prependTo=function(c){c=x(c);c.prepend(this);return this};c.before=function(c){if("set"==c.type){var a=this;c.forEach(function(c){var b=c.parent();a.node.parentNode.insertBefore(c.node,a.node);b&&b.add()});
this.parent().add();return this}c=x(c);var b=c.parent();this.node.parentNode.insertBefore(c.node,this.node);this.parent()&&this.parent().add();b&&b.add();c.paper=this.paper;return this};c.after=function(c){c=x(c);var a=c.parent();this.node.nextSibling?this.node.parentNode.insertBefore(c.node,this.node.nextSibling):this.node.parentNode.appendChild(c.node);this.parent()&&this.parent().add();a&&a.add();c.paper=this.paper;return this};c.insertBefore=function(c){c=x(c);var a=this.parent();c.node.parentNode.insertBefore(this.node,
c.node);this.paper=c.paper;a&&a.add();c.parent()&&c.parent().add();return this};c.insertAfter=function(c){c=x(c);var a=this.parent();c.node.parentNode.insertBefore(this.node,c.node.nextSibling);this.paper=c.paper;a&&a.add();c.parent()&&c.parent().add();return this};c.remove=function(){var c=this.parent();this.node.parentNode&&this.node.parentNode.removeChild(this.node);delete this.paper;this.removed=!0;c&&c.add();return this};c.select=function(c){return x(this.node.querySelector(c))};c.selectAll=
function(c){c=this.node.querySelectorAll(c);for(var b=(a.set||Array)(),m=0;m<c.length;m++)b.push(x(c[m]));return b};c.asPX=function(c,a){null==a&&(a=this.attr(c));return+b(this,c,a)};c.use=function(){var c,a=this.node.id;a||(a=this.id,v(this.node,{id:a}));c="linearGradient"==this.type||"radialGradient"==this.type||"pattern"==this.type?r(this.type,this.node.parentNode):r("use",this.node.parentNode);v(c.node,{"xlink:href":"#"+a});c.original=this;return c};var l=/\S+/g;c.addClass=function(c){var a=(c||
"").match(l)||[];c=this.node;var b=c.className.baseVal,m=b.match(l)||[],e,h,d;if(a.length){for(e=0;d=a[e++];)h=m.indexOf(d),~h||m.push(d);a=m.join(" ");b!=a&&(c.className.baseVal=a)}return this};c.removeClass=function(c){var a=(c||"").match(l)||[];c=this.node;var b=c.className.baseVal,m=b.match(l)||[],e,h;if(m.length){for(e=0;h=a[e++];)h=m.indexOf(h),~h&&m.splice(h,1);a=m.join(" ");b!=a&&(c.className.baseVal=a)}return this};c.hasClass=function(c){return!!~(this.node.className.baseVal.match(l)||[]).indexOf(c)};
c.toggleClass=function(c,a){if(null!=a)return a?this.addClass(c):this.removeClass(c);var b=(c||"").match(l)||[],m=this.node,e=m.className.baseVal,h=e.match(l)||[],d,f,E;for(d=0;E=b[d++];)f=h.indexOf(E),~f?h.splice(f,1):h.push(E);b=h.join(" ");e!=b&&(m.className.baseVal=b);return this};c.clone=function(){var c=x(this.node.cloneNode(!0));v(c.node,"id")&&v(c.node,{id:c.id});m(c);c.insertAfter(this);return c};c.toDefs=function(){u(this).appendChild(this.node);return this};c.pattern=c.toPattern=function(c,
a,b,m){var e=r("pattern",u(this));null==c&&(c=this.getBBox());y(c,"object")&&"x"in c&&(a=c.y,b=c.width,m=c.height,c=c.x);v(e.node,{x:c,y:a,width:b,height:m,patternUnits:"userSpaceOnUse",id:e.id,viewBox:[c,a,b,m].join(" ")});e.node.appendChild(this.node);return e};c.marker=function(c,a,b,m,e,h){var d=r("marker",u(this));null==c&&(c=this.getBBox());y(c,"object")&&"x"in c&&(a=c.y,b=c.width,m=c.height,e=c.refX||c.cx,h=c.refY||c.cy,c=c.x);v(d.node,{viewBox:[c,a,b,m].join(" "),markerWidth:b,markerHeight:m,
orient:"auto",refX:e||0,refY:h||0,id:d.id});d.node.appendChild(this.node);return d};var E=function(c,a,b,m){"function"!=typeof b||b.length||(m=b,b=L.linear);this.attr=c;this.dur=a;b&&(this.easing=b);m&&(this.callback=m)};a._.Animation=E;a.animation=function(c,a,b,m){return new E(c,a,b,m)};c.inAnim=function(){var c=[],a;for(a in this.anims)this.anims[h](a)&&function(a){c.push({anim:new E(a._attrs,a.dur,a.easing,a._callback),mina:a,curStatus:a.status(),status:function(c){return a.status(c)},stop:function(){a.stop()}})}(this.anims[a]);
return c};a.animate=function(c,a,b,m,e,h){"function"!=typeof e||e.length||(h=e,e=L.linear);var d=L.time();c=L(c,a,d,d+m,L.time,b,e);h&&k.once("mina.finish."+c.id,h);return c};c.stop=function(){for(var c=this.inAnim(),a=0,b=c.length;a<b;a++)c[a].stop();return this};c.animate=function(c,a,b,m){"function"!=typeof b||b.length||(m=b,b=L.linear);c instanceof E&&(m=c.callback,b=c.easing,a=b.dur,c=c.attr);var d=[],f=[],l={},t,ca,n,T=this,q;for(q in c)if(c[h](q)){T.equal?(n=T.equal(q,J(c[q])),t=n.from,ca=
n.to,n=n.f):(t=+T.attr(q),ca=+c[q]);var la=y(t,"array")?t.length:1;l[q]=e(d.length,d.length+la,n);d=d.concat(t);f=f.concat(ca)}t=L.time();var p=L(d,f,t,t+a,L.time,function(c){var a={},b;for(b in l)l[h](b)&&(a[b]=l[b](c));T.attr(a)},b);T.anims[p.id]=p;p._attrs=c;p._callback=m;k("snap.animcreated."+T.id,p);k.once("mina.finish."+p.id,function(){delete T.anims[p.id];m&&m.call(T)});k.once("mina.stop."+p.id,function(){delete T.anims[p.id]});return T};var T={};c.data=function(c,b){var m=T[this.id]=T[this.id]||
{};if(0==arguments.length)return k("snap.data.get."+this.id,this,m,null),m;if(1==arguments.length){if(a.is(c,"object")){for(var e in c)c[h](e)&&this.data(e,c[e]);return this}k("snap.data.get."+this.id,this,m[c],c);return m[c]}m[c]=b;k("snap.data.set."+this.id,this,b,c);return this};c.removeData=function(c){null==c?T[this.id]={}:T[this.id]&&delete T[this.id][c];return this};c.outerSVG=c.toString=d(1);c.innerSVG=d()})(e.prototype);a.parse=function(c){var a=G.doc.createDocumentFragment(),b=!0,m=G.doc.createElement("div");
c=J(c);c.match(/^\s*<\s*svg(?:\s|>)/)||(c="<svg>"+c+"</svg>",b=!1);m.innerHTML=c;if(c=m.getElementsByTagName("svg")[0])if(b)a=c;else for(;c.firstChild;)a.appendChild(c.firstChild);m.innerHTML=aa;return new l(a)};l.prototype.select=e.prototype.select;l.prototype.selectAll=e.prototype.selectAll;a.fragment=function(){for(var c=Array.prototype.slice.call(arguments,0),b=G.doc.createDocumentFragment(),m=0,e=c.length;m<e;m++){var h=c[m];h.node&&h.node.nodeType&&b.appendChild(h.node);h.nodeType&&b.appendChild(h);
"string"==typeof h&&b.appendChild(a.parse(h).node)}return new l(b)};a._.make=r;a._.wrap=x;s.prototype.el=function(c,a){var b=r(c,this.node);a&&b.attr(a);return b};k.on("snap.util.getattr",function(){var c=k.nt(),c=c.substring(c.lastIndexOf(".")+1),a=c.replace(/[A-Z]/g,function(c){return"-"+c.toLowerCase()});return pa[h](a)?this.node.ownerDocument.defaultView.getComputedStyle(this.node,null).getPropertyValue(a):v(this.node,c)});var pa={"alignment-baseline":0,"baseline-shift":0,clip:0,"clip-path":0,
"clip-rule":0,color:0,"color-interpolation":0,"color-interpolation-filters":0,"color-profile":0,"color-rendering":0,cursor:0,direction:0,display:0,"dominant-baseline":0,"enable-background":0,fill:0,"fill-opacity":0,"fill-rule":0,filter:0,"flood-color":0,"flood-opacity":0,font:0,"font-family":0,"font-size":0,"font-size-adjust":0,"font-stretch":0,"font-style":0,"font-variant":0,"font-weight":0,"glyph-orientation-horizontal":0,"glyph-orientation-vertical":0,"image-rendering":0,kerning:0,"letter-spacing":0,
"lighting-color":0,marker:0,"marker-end":0,"marker-mid":0,"marker-start":0,mask:0,opacity:0,overflow:0,"pointer-events":0,"shape-rendering":0,"stop-color":0,"stop-opacity":0,stroke:0,"stroke-dasharray":0,"stroke-dashoffset":0,"stroke-linecap":0,"stroke-linejoin":0,"stroke-miterlimit":0,"stroke-opacity":0,"stroke-width":0,"text-anchor":0,"text-decoration":0,"text-rendering":0,"unicode-bidi":0,visibility:0,"word-spacing":0,"writing-mode":0};k.on("snap.util.attr",function(c){var a=k.nt(),b={},a=a.substring(a.lastIndexOf(".")+
1);b[a]=c;var m=a.replace(/-(\w)/gi,function(c,a){return a.toUpperCase()}),a=a.replace(/[A-Z]/g,function(c){return"-"+c.toLowerCase()});pa[h](a)?this.node.style[m]=null==c?aa:c:v(this.node,b)});a.ajax=function(c,a,b,m){var e=new XMLHttpRequest,h=V();if(e){if(y(a,"function"))m=b,b=a,a=null;else if(y(a,"object")){var d=[],f;for(f in a)a.hasOwnProperty(f)&&d.push(encodeURIComponent(f)+"="+encodeURIComponent(a[f]));a=d.join("&")}e.open(a?"POST":"GET",c,!0);a&&(e.setRequestHeader("X-Requested-With","XMLHttpRequest"),
e.setRequestHeader("Content-type","application/x-www-form-urlencoded"));b&&(k.once("snap.ajax."+h+".0",b),k.once("snap.ajax."+h+".200",b),k.once("snap.ajax."+h+".304",b));e.onreadystatechange=function(){4==e.readyState&&k("snap.ajax."+h+"."+e.status,m,e)};if(4==e.readyState)return e;e.send(a);return e}};a.load=function(c,b,m){a.ajax(c,function(c){c=a.parse(c.responseText);m?b.call(m,c):b(c)})};a.getElementByPoint=function(c,a){var b,m,e=G.doc.elementFromPoint(c,a);if(G.win.opera&&"svg"==e.tagName){b=
e;m=b.getBoundingClientRect();b=b.ownerDocument;var h=b.body,d=b.documentElement;b=m.top+(g.win.pageYOffset||d.scrollTop||h.scrollTop)-(d.clientTop||h.clientTop||0);m=m.left+(g.win.pageXOffset||d.scrollLeft||h.scrollLeft)-(d.clientLeft||h.clientLeft||0);h=e.createSVGRect();h.x=c-m;h.y=a-b;h.width=h.height=1;b=e.getIntersectionList(h,null);b.length&&(e=b[b.length-1])}return e?x(e):null};a.plugin=function(c){c(a,e,s,G,l)};return G.win.Snap=a}();C.plugin(function(a,k,y,M,A){function w(a,d,f,b,q,e){null==
d&&"[object SVGMatrix]"==z.call(a)?(this.a=a.a,this.b=a.b,this.c=a.c,this.d=a.d,this.e=a.e,this.f=a.f):null!=a?(this.a=+a,this.b=+d,this.c=+f,this.d=+b,this.e=+q,this.f=+e):(this.a=1,this.c=this.b=0,this.d=1,this.f=this.e=0)}var z=Object.prototype.toString,d=String,f=Math;(function(n){function k(a){return a[0]*a[0]+a[1]*a[1]}function p(a){var d=f.sqrt(k(a));a[0]&&(a[0]/=d);a[1]&&(a[1]/=d)}n.add=function(a,d,e,f,n,p){var k=[[],[],[] ],u=[[this.a,this.c,this.e],[this.b,this.d,this.f],[0,0,1] ];d=[[a,
e,n],[d,f,p],[0,0,1] ];a&&a instanceof w&&(d=[[a.a,a.c,a.e],[a.b,a.d,a.f],[0,0,1] ]);for(a=0;3>a;a++)for(e=0;3>e;e++){for(f=n=0;3>f;f++)n+=u[a][f]*d[f][e];k[a][e]=n}this.a=k[0][0];this.b=k[1][0];this.c=k[0][1];this.d=k[1][1];this.e=k[0][2];this.f=k[1][2];return this};n.invert=function(){var a=this.a*this.d-this.b*this.c;return new w(this.d/a,-this.b/a,-this.c/a,this.a/a,(this.c*this.f-this.d*this.e)/a,(this.b*this.e-this.a*this.f)/a)};n.clone=function(){return new w(this.a,this.b,this.c,this.d,this.e,
this.f)};n.translate=function(a,d){return this.add(1,0,0,1,a,d)};n.scale=function(a,d,e,f){null==d&&(d=a);(e||f)&&this.add(1,0,0,1,e,f);this.add(a,0,0,d,0,0);(e||f)&&this.add(1,0,0,1,-e,-f);return this};n.rotate=function(b,d,e){b=a.rad(b);d=d||0;e=e||0;var l=+f.cos(b).toFixed(9);b=+f.sin(b).toFixed(9);this.add(l,b,-b,l,d,e);return this.add(1,0,0,1,-d,-e)};n.x=function(a,d){return a*this.a+d*this.c+this.e};n.y=function(a,d){return a*this.b+d*this.d+this.f};n.get=function(a){return+this[d.fromCharCode(97+
a)].toFixed(4)};n.toString=function(){return"matrix("+[this.get(0),this.get(1),this.get(2),this.get(3),this.get(4),this.get(5)].join()+")"};n.offset=function(){return[this.e.toFixed(4),this.f.toFixed(4)]};n.determinant=function(){return this.a*this.d-this.b*this.c};n.split=function(){var b={};b.dx=this.e;b.dy=this.f;var d=[[this.a,this.c],[this.b,this.d] ];b.scalex=f.sqrt(k(d[0]));p(d[0]);b.shear=d[0][0]*d[1][0]+d[0][1]*d[1][1];d[1]=[d[1][0]-d[0][0]*b.shear,d[1][1]-d[0][1]*b.shear];b.scaley=f.sqrt(k(d[1]));
p(d[1]);b.shear/=b.scaley;0>this.determinant()&&(b.scalex=-b.scalex);var e=-d[0][1],d=d[1][1];0>d?(b.rotate=a.deg(f.acos(d)),0>e&&(b.rotate=360-b.rotate)):b.rotate=a.deg(f.asin(e));b.isSimple=!+b.shear.toFixed(9)&&(b.scalex.toFixed(9)==b.scaley.toFixed(9)||!b.rotate);b.isSuperSimple=!+b.shear.toFixed(9)&&b.scalex.toFixed(9)==b.scaley.toFixed(9)&&!b.rotate;b.noRotation=!+b.shear.toFixed(9)&&!b.rotate;return b};n.toTransformString=function(a){a=a||this.split();if(+a.shear.toFixed(9))return"m"+[this.get(0),
this.get(1),this.get(2),this.get(3),this.get(4),this.get(5)];a.scalex=+a.scalex.toFixed(4);a.scaley=+a.scaley.toFixed(4);a.rotate=+a.rotate.toFixed(4);return(a.dx||a.dy?"t"+[+a.dx.toFixed(4),+a.dy.toFixed(4)]:"")+(1!=a.scalex||1!=a.scaley?"s"+[a.scalex,a.scaley,0,0]:"")+(a.rotate?"r"+[+a.rotate.toFixed(4),0,0]:"")}})(w.prototype);a.Matrix=w;a.matrix=function(a,d,f,b,k,e){return new w(a,d,f,b,k,e)}});C.plugin(function(a,v,y,M,A){function w(h){return function(d){k.stop();d instanceof A&&1==d.node.childNodes.length&&
("radialGradient"==d.node.firstChild.tagName||"linearGradient"==d.node.firstChild.tagName||"pattern"==d.node.firstChild.tagName)&&(d=d.node.firstChild,b(this).appendChild(d),d=u(d));if(d instanceof v)if("radialGradient"==d.type||"linearGradient"==d.type||"pattern"==d.type){d.node.id||e(d.node,{id:d.id});var f=l(d.node.id)}else f=d.attr(h);else f=a.color(d),f.error?(f=a(b(this).ownerSVGElement).gradient(d))?(f.node.id||e(f.node,{id:f.id}),f=l(f.node.id)):f=d:f=r(f);d={};d[h]=f;e(this.node,d);this.node.style[h]=
x}}function z(a){k.stop();a==+a&&(a+="px");this.node.style.fontSize=a}function d(a){var b=[];a=a.childNodes;for(var e=0,f=a.length;e<f;e++){var l=a[e];3==l.nodeType&&b.push(l.nodeValue);"tspan"==l.tagName&&(1==l.childNodes.length&&3==l.firstChild.nodeType?b.push(l.firstChild.nodeValue):b.push(d(l)))}return b}function f(){k.stop();return this.node.style.fontSize}var n=a._.make,u=a._.wrap,p=a.is,b=a._.getSomeDefs,q=/^url\(#?([^)]+)\)$/,e=a._.$,l=a.url,r=String,s=a._.separator,x="";k.on("snap.util.attr.mask",
function(a){if(a instanceof v||a instanceof A){k.stop();a instanceof A&&1==a.node.childNodes.length&&(a=a.node.firstChild,b(this).appendChild(a),a=u(a));if("mask"==a.type)var d=a;else d=n("mask",b(this)),d.node.appendChild(a.node);!d.node.id&&e(d.node,{id:d.id});e(this.node,{mask:l(d.id)})}});(function(a){k.on("snap.util.attr.clip",a);k.on("snap.util.attr.clip-path",a);k.on("snap.util.attr.clipPath",a)})(function(a){if(a instanceof v||a instanceof A){k.stop();if("clipPath"==a.type)var d=a;else d=
n("clipPath",b(this)),d.node.appendChild(a.node),!d.node.id&&e(d.node,{id:d.id});e(this.node,{"clip-path":l(d.id)})}});k.on("snap.util.attr.fill",w("fill"));k.on("snap.util.attr.stroke",w("stroke"));var G=/^([lr])(?:\(([^)]*)\))?(.*)$/i;k.on("snap.util.grad.parse",function(a){a=r(a);var b=a.match(G);if(!b)return null;a=b[1];var e=b[2],b=b[3],e=e.split(/\s*,\s*/).map(function(a){return+a==a?+a:a});1==e.length&&0==e[0]&&(e=[]);b=b.split("-");b=b.map(function(a){a=a.split(":");var b={color:a[0]};a[1]&&
(b.offset=parseFloat(a[1]));return b});return{type:a,params:e,stops:b}});k.on("snap.util.attr.d",function(b){k.stop();p(b,"array")&&p(b[0],"array")&&(b=a.path.toString.call(b));b=r(b);b.match(/[ruo]/i)&&(b=a.path.toAbsolute(b));e(this.node,{d:b})})(-1);k.on("snap.util.attr.#text",function(a){k.stop();a=r(a);for(a=M.doc.createTextNode(a);this.node.firstChild;)this.node.removeChild(this.node.firstChild);this.node.appendChild(a)})(-1);k.on("snap.util.attr.path",function(a){k.stop();this.attr({d:a})})(-1);
k.on("snap.util.attr.class",function(a){k.stop();this.node.className.baseVal=a})(-1);k.on("snap.util.attr.viewBox",function(a){a=p(a,"object")&&"x"in a?[a.x,a.y,a.width,a.height].join(" "):p(a,"array")?a.join(" "):a;e(this.node,{viewBox:a});k.stop()})(-1);k.on("snap.util.attr.transform",function(a){this.transform(a);k.stop()})(-1);k.on("snap.util.attr.r",function(a){"rect"==this.type&&(k.stop(),e(this.node,{rx:a,ry:a}))})(-1);k.on("snap.util.attr.textpath",function(a){k.stop();if("text"==this.type){var d,
f;if(!a&&this.textPath){for(a=this.textPath;a.node.firstChild;)this.node.appendChild(a.node.firstChild);a.remove();delete this.textPath}else if(p(a,"string")?(d=b(this),a=u(d.parentNode).path(a),d.appendChild(a.node),d=a.id,a.attr({id:d})):(a=u(a),a instanceof v&&(d=a.attr("id"),d||(d=a.id,a.attr({id:d})))),d)if(a=this.textPath,f=this.node,a)a.attr({"xlink:href":"#"+d});else{for(a=e("textPath",{"xlink:href":"#"+d});f.firstChild;)a.appendChild(f.firstChild);f.appendChild(a);this.textPath=u(a)}}})(-1);
k.on("snap.util.attr.text",function(a){if("text"==this.type){for(var b=this.node,d=function(a){var b=e("tspan");if(p(a,"array"))for(var f=0;f<a.length;f++)b.appendChild(d(a[f]));else b.appendChild(M.doc.createTextNode(a));b.normalize&&b.normalize();return b};b.firstChild;)b.removeChild(b.firstChild);for(a=d(a);a.firstChild;)b.appendChild(a.firstChild)}k.stop()})(-1);k.on("snap.util.attr.fontSize",z)(-1);k.on("snap.util.attr.font-size",z)(-1);k.on("snap.util.getattr.transform",function(){k.stop();
return this.transform()})(-1);k.on("snap.util.getattr.textpath",function(){k.stop();return this.textPath})(-1);(function(){function b(d){return function(){k.stop();var b=M.doc.defaultView.getComputedStyle(this.node,null).getPropertyValue("marker-"+d);return"none"==b?b:a(M.doc.getElementById(b.match(q)[1]))}}function d(a){return function(b){k.stop();var d="marker"+a.charAt(0).toUpperCase()+a.substring(1);if(""==b||!b)this.node.style[d]="none";else if("marker"==b.type){var f=b.node.id;f||e(b.node,{id:b.id});
this.node.style[d]=l(f)}}}k.on("snap.util.getattr.marker-end",b("end"))(-1);k.on("snap.util.getattr.markerEnd",b("end"))(-1);k.on("snap.util.getattr.marker-start",b("start"))(-1);k.on("snap.util.getattr.markerStart",b("start"))(-1);k.on("snap.util.getattr.marker-mid",b("mid"))(-1);k.on("snap.util.getattr.markerMid",b("mid"))(-1);k.on("snap.util.attr.marker-end",d("end"))(-1);k.on("snap.util.attr.markerEnd",d("end"))(-1);k.on("snap.util.attr.marker-start",d("start"))(-1);k.on("snap.util.attr.markerStart",
d("start"))(-1);k.on("snap.util.attr.marker-mid",d("mid"))(-1);k.on("snap.util.attr.markerMid",d("mid"))(-1)})();k.on("snap.util.getattr.r",function(){if("rect"==this.type&&e(this.node,"rx")==e(this.node,"ry"))return k.stop(),e(this.node,"rx")})(-1);k.on("snap.util.getattr.text",function(){if("text"==this.type||"tspan"==this.type){k.stop();var a=d(this.node);return 1==a.length?a[0]:a}})(-1);k.on("snap.util.getattr.#text",function(){return this.node.textContent})(-1);k.on("snap.util.getattr.viewBox",
function(){k.stop();var b=e(this.node,"viewBox");if(b)return b=b.split(s),a._.box(+b[0],+b[1],+b[2],+b[3])})(-1);k.on("snap.util.getattr.points",function(){var a=e(this.node,"points");k.stop();if(a)return a.split(s)})(-1);k.on("snap.util.getattr.path",function(){var a=e(this.node,"d");k.stop();return a})(-1);k.on("snap.util.getattr.class",function(){return this.node.className.baseVal})(-1);k.on("snap.util.getattr.fontSize",f)(-1);k.on("snap.util.getattr.font-size",f)(-1)});C.plugin(function(a,v,y,
M,A){function w(a){return a}function z(a){return function(b){return+b.toFixed(3)+a}}var d={"+":function(a,b){return a+b},"-":function(a,b){return a-b},"/":function(a,b){return a/b},"*":function(a,b){return a*b}},f=String,n=/[a-z]+$/i,u=/^\s*([+\-\/*])\s*=\s*([\d.eE+\-]+)\s*([^\d\s]+)?\s*$/;k.on("snap.util.attr",function(a){if(a=f(a).match(u)){var b=k.nt(),b=b.substring(b.lastIndexOf(".")+1),q=this.attr(b),e={};k.stop();var l=a[3]||"",r=q.match(n),s=d[a[1] ];r&&r==l?a=s(parseFloat(q),+a[2]):(q=this.asPX(b),
a=s(this.asPX(b),this.asPX(b,a[2]+l)));isNaN(q)||isNaN(a)||(e[b]=a,this.attr(e))}})(-10);k.on("snap.util.equal",function(a,b){var q=f(this.attr(a)||""),e=f(b).match(u);if(e){k.stop();var l=e[3]||"",r=q.match(n),s=d[e[1] ];if(r&&r==l)return{from:parseFloat(q),to:s(parseFloat(q),+e[2]),f:z(r)};q=this.asPX(a);return{from:q,to:s(q,this.asPX(a,e[2]+l)),f:w}}})(-10)});C.plugin(function(a,v,y,M,A){var w=y.prototype,z=a.is;w.rect=function(a,d,k,p,b,q){var e;null==q&&(q=b);z(a,"object")&&"[object Object]"==
a?e=a:null!=a&&(e={x:a,y:d,width:k,height:p},null!=b&&(e.rx=b,e.ry=q));return this.el("rect",e)};w.circle=function(a,d,k){var p;z(a,"object")&&"[object Object]"==a?p=a:null!=a&&(p={cx:a,cy:d,r:k});return this.el("circle",p)};var d=function(){function a(){this.parentNode.removeChild(this)}return function(d,k){var p=M.doc.createElement("img"),b=M.doc.body;p.style.cssText="position:absolute;left:-9999em;top:-9999em";p.onload=function(){k.call(p);p.onload=p.onerror=null;b.removeChild(p)};p.onerror=a;
b.appendChild(p);p.src=d}}();w.image=function(f,n,k,p,b){var q=this.el("image");if(z(f,"object")&&"src"in f)q.attr(f);else if(null!=f){var e={"xlink:href":f,preserveAspectRatio:"none"};null!=n&&null!=k&&(e.x=n,e.y=k);null!=p&&null!=b?(e.width=p,e.height=b):d(f,function(){a._.$(q.node,{width:this.offsetWidth,height:this.offsetHeight})});a._.$(q.node,e)}return q};w.ellipse=function(a,d,k,p){var b;z(a,"object")&&"[object Object]"==a?b=a:null!=a&&(b={cx:a,cy:d,rx:k,ry:p});return this.el("ellipse",b)};
w.path=function(a){var d;z(a,"object")&&!z(a,"array")?d=a:a&&(d={d:a});return this.el("path",d)};w.group=w.g=function(a){var d=this.el("g");1==arguments.length&&a&&!a.type?d.attr(a):arguments.length&&d.add(Array.prototype.slice.call(arguments,0));return d};w.svg=function(a,d,k,p,b,q,e,l){var r={};z(a,"object")&&null==d?r=a:(null!=a&&(r.x=a),null!=d&&(r.y=d),null!=k&&(r.width=k),null!=p&&(r.height=p),null!=b&&null!=q&&null!=e&&null!=l&&(r.viewBox=[b,q,e,l]));return this.el("svg",r)};w.mask=function(a){var d=
this.el("mask");1==arguments.length&&a&&!a.type?d.attr(a):arguments.length&&d.add(Array.prototype.slice.call(arguments,0));return d};w.ptrn=function(a,d,k,p,b,q,e,l){if(z(a,"object"))var r=a;else arguments.length?(r={},null!=a&&(r.x=a),null!=d&&(r.y=d),null!=k&&(r.width=k),null!=p&&(r.height=p),null!=b&&null!=q&&null!=e&&null!=l&&(r.viewBox=[b,q,e,l])):r={patternUnits:"userSpaceOnUse"};return this.el("pattern",r)};w.use=function(a){return null!=a?(make("use",this.node),a instanceof v&&(a.attr("id")||
a.attr({id:ID()}),a=a.attr("id")),this.el("use",{"xlink:href":a})):v.prototype.use.call(this)};w.text=function(a,d,k){var p={};z(a,"object")?p=a:null!=a&&(p={x:a,y:d,text:k||""});return this.el("text",p)};w.line=function(a,d,k,p){var b={};z(a,"object")?b=a:null!=a&&(b={x1:a,x2:k,y1:d,y2:p});return this.el("line",b)};w.polyline=function(a){1<arguments.length&&(a=Array.prototype.slice.call(arguments,0));var d={};z(a,"object")&&!z(a,"array")?d=a:null!=a&&(d={points:a});return this.el("polyline",d)};
w.polygon=function(a){1<arguments.length&&(a=Array.prototype.slice.call(arguments,0));var d={};z(a,"object")&&!z(a,"array")?d=a:null!=a&&(d={points:a});return this.el("polygon",d)};(function(){function d(){return this.selectAll("stop")}function n(b,d){var f=e("stop"),k={offset:+d+"%"};b=a.color(b);k["stop-color"]=b.hex;1>b.opacity&&(k["stop-opacity"]=b.opacity);e(f,k);this.node.appendChild(f);return this}function u(){if("linearGradient"==this.type){var b=e(this.node,"x1")||0,d=e(this.node,"x2")||
1,f=e(this.node,"y1")||0,k=e(this.node,"y2")||0;return a._.box(b,f,math.abs(d-b),math.abs(k-f))}b=this.node.r||0;return a._.box((this.node.cx||0.5)-b,(this.node.cy||0.5)-b,2*b,2*b)}function p(a,d){function f(a,b){for(var d=(b-u)/(a-w),e=w;e<a;e++)h[e].offset=+(+u+d*(e-w)).toFixed(2);w=a;u=b}var n=k("snap.util.grad.parse",null,d).firstDefined(),p;if(!n)return null;n.params.unshift(a);p="l"==n.type.toLowerCase()?b.apply(0,n.params):q.apply(0,n.params);n.type!=n.type.toLowerCase()&&e(p.node,{gradientUnits:"userSpaceOnUse"});
var h=n.stops,n=h.length,u=0,w=0;n--;for(var v=0;v<n;v++)"offset"in h[v]&&f(v,h[v].offset);h[n].offset=h[n].offset||100;f(n,h[n].offset);for(v=0;v<=n;v++){var y=h[v];p.addStop(y.color,y.offset)}return p}function b(b,k,p,q,w){b=a._.make("linearGradient",b);b.stops=d;b.addStop=n;b.getBBox=u;null!=k&&e(b.node,{x1:k,y1:p,x2:q,y2:w});return b}function q(b,k,p,q,w,h){b=a._.make("radialGradient",b);b.stops=d;b.addStop=n;b.getBBox=u;null!=k&&e(b.node,{cx:k,cy:p,r:q});null!=w&&null!=h&&e(b.node,{fx:w,fy:h});
return b}var e=a._.$;w.gradient=function(a){return p(this.defs,a)};w.gradientLinear=function(a,d,e,f){return b(this.defs,a,d,e,f)};w.gradientRadial=function(a,b,d,e,f){return q(this.defs,a,b,d,e,f)};w.toString=function(){var b=this.node.ownerDocument,d=b.createDocumentFragment(),b=b.createElement("div"),e=this.node.cloneNode(!0);d.appendChild(b);b.appendChild(e);a._.$(e,{xmlns:"http://www.w3.org/2000/svg"});b=b.innerHTML;d.removeChild(d.firstChild);return b};w.clear=function(){for(var a=this.node.firstChild,
b;a;)b=a.nextSibling,"defs"!=a.tagName?a.parentNode.removeChild(a):w.clear.call({node:a}),a=b}})()});C.plugin(function(a,k,y,M){function A(a){var b=A.ps=A.ps||{};b[a]?b[a].sleep=100:b[a]={sleep:100};setTimeout(function(){for(var d in b)b[L](d)&&d!=a&&(b[d].sleep--,!b[d].sleep&&delete b[d])});return b[a]}function w(a,b,d,e){null==a&&(a=b=d=e=0);null==b&&(b=a.y,d=a.width,e=a.height,a=a.x);return{x:a,y:b,width:d,w:d,height:e,h:e,x2:a+d,y2:b+e,cx:a+d/2,cy:b+e/2,r1:F.min(d,e)/2,r2:F.max(d,e)/2,r0:F.sqrt(d*
d+e*e)/2,path:s(a,b,d,e),vb:[a,b,d,e].join(" ")}}function z(){return this.join(",").replace(N,"$1")}function d(a){a=C(a);a.toString=z;return a}function f(a,b,d,h,f,k,l,n,p){if(null==p)return e(a,b,d,h,f,k,l,n);if(0>p||e(a,b,d,h,f,k,l,n)<p)p=void 0;else{var q=0.5,O=1-q,s;for(s=e(a,b,d,h,f,k,l,n,O);0.01<Z(s-p);)q/=2,O+=(s<p?1:-1)*q,s=e(a,b,d,h,f,k,l,n,O);p=O}return u(a,b,d,h,f,k,l,n,p)}function n(b,d){function e(a){return+(+a).toFixed(3)}return a._.cacher(function(a,h,l){a instanceof k&&(a=a.attr("d"));
a=I(a);for(var n,p,D,q,O="",s={},c=0,t=0,r=a.length;t<r;t++){D=a[t];if("M"==D[0])n=+D[1],p=+D[2];else{q=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6]);if(c+q>h){if(d&&!s.start){n=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6],h-c);O+=["C"+e(n.start.x),e(n.start.y),e(n.m.x),e(n.m.y),e(n.x),e(n.y)];if(l)return O;s.start=O;O=["M"+e(n.x),e(n.y)+"C"+e(n.n.x),e(n.n.y),e(n.end.x),e(n.end.y),e(D[5]),e(D[6])].join();c+=q;n=+D[5];p=+D[6];continue}if(!b&&!d)return n=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6],h-c)}c+=q;n=+D[5];p=+D[6]}O+=
D.shift()+D}s.end=O;return n=b?c:d?s:u(n,p,D[0],D[1],D[2],D[3],D[4],D[5],1)},null,a._.clone)}function u(a,b,d,e,h,f,k,l,n){var p=1-n,q=ma(p,3),s=ma(p,2),c=n*n,t=c*n,r=q*a+3*s*n*d+3*p*n*n*h+t*k,q=q*b+3*s*n*e+3*p*n*n*f+t*l,s=a+2*n*(d-a)+c*(h-2*d+a),t=b+2*n*(e-b)+c*(f-2*e+b),x=d+2*n*(h-d)+c*(k-2*h+d),c=e+2*n*(f-e)+c*(l-2*f+e);a=p*a+n*d;b=p*b+n*e;h=p*h+n*k;f=p*f+n*l;l=90-180*F.atan2(s-x,t-c)/S;return{x:r,y:q,m:{x:s,y:t},n:{x:x,y:c},start:{x:a,y:b},end:{x:h,y:f},alpha:l}}function p(b,d,e,h,f,n,k,l){a.is(b,
"array")||(b=[b,d,e,h,f,n,k,l]);b=U.apply(null,b);return w(b.min.x,b.min.y,b.max.x-b.min.x,b.max.y-b.min.y)}function b(a,b,d){return b>=a.x&&b<=a.x+a.width&&d>=a.y&&d<=a.y+a.height}function q(a,d){a=w(a);d=w(d);return b(d,a.x,a.y)||b(d,a.x2,a.y)||b(d,a.x,a.y2)||b(d,a.x2,a.y2)||b(a,d.x,d.y)||b(a,d.x2,d.y)||b(a,d.x,d.y2)||b(a,d.x2,d.y2)||(a.x<d.x2&&a.x>d.x||d.x<a.x2&&d.x>a.x)&&(a.y<d.y2&&a.y>d.y||d.y<a.y2&&d.y>a.y)}function e(a,b,d,e,h,f,n,k,l){null==l&&(l=1);l=(1<l?1:0>l?0:l)/2;for(var p=[-0.1252,
0.1252,-0.3678,0.3678,-0.5873,0.5873,-0.7699,0.7699,-0.9041,0.9041,-0.9816,0.9816],q=[0.2491,0.2491,0.2335,0.2335,0.2032,0.2032,0.1601,0.1601,0.1069,0.1069,0.0472,0.0472],s=0,c=0;12>c;c++)var t=l*p[c]+l,r=t*(t*(-3*a+9*d-9*h+3*n)+6*a-12*d+6*h)-3*a+3*d,t=t*(t*(-3*b+9*e-9*f+3*k)+6*b-12*e+6*f)-3*b+3*e,s=s+q[c]*F.sqrt(r*r+t*t);return l*s}function l(a,b,d){a=I(a);b=I(b);for(var h,f,l,n,k,s,r,O,x,c,t=d?0:[],w=0,v=a.length;w<v;w++)if(x=a[w],"M"==x[0])h=k=x[1],f=s=x[2];else{"C"==x[0]?(x=[h,f].concat(x.slice(1)),
h=x[6],f=x[7]):(x=[h,f,h,f,k,s,k,s],h=k,f=s);for(var G=0,y=b.length;G<y;G++)if(c=b[G],"M"==c[0])l=r=c[1],n=O=c[2];else{"C"==c[0]?(c=[l,n].concat(c.slice(1)),l=c[6],n=c[7]):(c=[l,n,l,n,r,O,r,O],l=r,n=O);var z;var K=x,B=c;z=d;var H=p(K),J=p(B);if(q(H,J)){for(var H=e.apply(0,K),J=e.apply(0,B),H=~~(H/8),J=~~(J/8),U=[],A=[],F={},M=z?0:[],P=0;P<H+1;P++){var C=u.apply(0,K.concat(P/H));U.push({x:C.x,y:C.y,t:P/H})}for(P=0;P<J+1;P++)C=u.apply(0,B.concat(P/J)),A.push({x:C.x,y:C.y,t:P/J});for(P=0;P<H;P++)for(K=
0;K<J;K++){var Q=U[P],L=U[P+1],B=A[K],C=A[K+1],N=0.001>Z(L.x-Q.x)?"y":"x",S=0.001>Z(C.x-B.x)?"y":"x",R;R=Q.x;var Y=Q.y,V=L.x,ea=L.y,fa=B.x,ga=B.y,ha=C.x,ia=C.y;if(W(R,V)<X(fa,ha)||X(R,V)>W(fa,ha)||W(Y,ea)<X(ga,ia)||X(Y,ea)>W(ga,ia))R=void 0;else{var $=(R*ea-Y*V)*(fa-ha)-(R-V)*(fa*ia-ga*ha),aa=(R*ea-Y*V)*(ga-ia)-(Y-ea)*(fa*ia-ga*ha),ja=(R-V)*(ga-ia)-(Y-ea)*(fa-ha);if(ja){var $=$/ja,aa=aa/ja,ja=+$.toFixed(2),ba=+aa.toFixed(2);R=ja<+X(R,V).toFixed(2)||ja>+W(R,V).toFixed(2)||ja<+X(fa,ha).toFixed(2)||
ja>+W(fa,ha).toFixed(2)||ba<+X(Y,ea).toFixed(2)||ba>+W(Y,ea).toFixed(2)||ba<+X(ga,ia).toFixed(2)||ba>+W(ga,ia).toFixed(2)?void 0:{x:$,y:aa}}else R=void 0}R&&F[R.x.toFixed(4)]!=R.y.toFixed(4)&&(F[R.x.toFixed(4)]=R.y.toFixed(4),Q=Q.t+Z((R[N]-Q[N])/(L[N]-Q[N]))*(L.t-Q.t),B=B.t+Z((R[S]-B[S])/(C[S]-B[S]))*(C.t-B.t),0<=Q&&1>=Q&&0<=B&&1>=B&&(z?M++:M.push({x:R.x,y:R.y,t1:Q,t2:B})))}z=M}else z=z?0:[];if(d)t+=z;else{H=0;for(J=z.length;H<J;H++)z[H].segment1=w,z[H].segment2=G,z[H].bez1=x,z[H].bez2=c;t=t.concat(z)}}}return t}
function r(a){var b=A(a);if(b.bbox)return C(b.bbox);if(!a)return w();a=I(a);for(var d=0,e=0,h=[],f=[],l,n=0,k=a.length;n<k;n++)l=a[n],"M"==l[0]?(d=l[1],e=l[2],h.push(d),f.push(e)):(d=U(d,e,l[1],l[2],l[3],l[4],l[5],l[6]),h=h.concat(d.min.x,d.max.x),f=f.concat(d.min.y,d.max.y),d=l[5],e=l[6]);a=X.apply(0,h);l=X.apply(0,f);h=W.apply(0,h);f=W.apply(0,f);f=w(a,l,h-a,f-l);b.bbox=C(f);return f}function s(a,b,d,e,h){if(h)return[["M",+a+ +h,b],["l",d-2*h,0],["a",h,h,0,0,1,h,h],["l",0,e-2*h],["a",h,h,0,0,1,
-h,h],["l",2*h-d,0],["a",h,h,0,0,1,-h,-h],["l",0,2*h-e],["a",h,h,0,0,1,h,-h],["z"] ];a=[["M",a,b],["l",d,0],["l",0,e],["l",-d,0],["z"] ];a.toString=z;return a}function x(a,b,d,e,h){null==h&&null==e&&(e=d);a=+a;b=+b;d=+d;e=+e;if(null!=h){var f=Math.PI/180,l=a+d*Math.cos(-e*f);a+=d*Math.cos(-h*f);var n=b+d*Math.sin(-e*f);b+=d*Math.sin(-h*f);d=[["M",l,n],["A",d,d,0,+(180<h-e),0,a,b] ]}else d=[["M",a,b],["m",0,-e],["a",d,e,0,1,1,0,2*e],["a",d,e,0,1,1,0,-2*e],["z"] ];d.toString=z;return d}function G(b){var e=
A(b);if(e.abs)return d(e.abs);Q(b,"array")&&Q(b&&b[0],"array")||(b=a.parsePathString(b));if(!b||!b.length)return[["M",0,0] ];var h=[],f=0,l=0,n=0,k=0,p=0;"M"==b[0][0]&&(f=+b[0][1],l=+b[0][2],n=f,k=l,p++,h[0]=["M",f,l]);for(var q=3==b.length&&"M"==b[0][0]&&"R"==b[1][0].toUpperCase()&&"Z"==b[2][0].toUpperCase(),s,r,w=p,c=b.length;w<c;w++){h.push(s=[]);r=b[w];p=r[0];if(p!=p.toUpperCase())switch(s[0]=p.toUpperCase(),s[0]){case "A":s[1]=r[1];s[2]=r[2];s[3]=r[3];s[4]=r[4];s[5]=r[5];s[6]=+r[6]+f;s[7]=+r[7]+
l;break;case "V":s[1]=+r[1]+l;break;case "H":s[1]=+r[1]+f;break;case "R":for(var t=[f,l].concat(r.slice(1)),u=2,v=t.length;u<v;u++)t[u]=+t[u]+f,t[++u]=+t[u]+l;h.pop();h=h.concat(P(t,q));break;case "O":h.pop();t=x(f,l,r[1],r[2]);t.push(t[0]);h=h.concat(t);break;case "U":h.pop();h=h.concat(x(f,l,r[1],r[2],r[3]));s=["U"].concat(h[h.length-1].slice(-2));break;case "M":n=+r[1]+f,k=+r[2]+l;default:for(u=1,v=r.length;u<v;u++)s[u]=+r[u]+(u%2?f:l)}else if("R"==p)t=[f,l].concat(r.slice(1)),h.pop(),h=h.concat(P(t,
q)),s=["R"].concat(r.slice(-2));else if("O"==p)h.pop(),t=x(f,l,r[1],r[2]),t.push(t[0]),h=h.concat(t);else if("U"==p)h.pop(),h=h.concat(x(f,l,r[1],r[2],r[3])),s=["U"].concat(h[h.length-1].slice(-2));else for(t=0,u=r.length;t<u;t++)s[t]=r[t];p=p.toUpperCase();if("O"!=p)switch(s[0]){case "Z":f=+n;l=+k;break;case "H":f=s[1];break;case "V":l=s[1];break;case "M":n=s[s.length-2],k=s[s.length-1];default:f=s[s.length-2],l=s[s.length-1]}}h.toString=z;e.abs=d(h);return h}function h(a,b,d,e){return[a,b,d,e,d,
e]}function J(a,b,d,e,h,f){var l=1/3,n=2/3;return[l*a+n*d,l*b+n*e,l*h+n*d,l*f+n*e,h,f]}function K(b,d,e,h,f,l,n,k,p,s){var r=120*S/180,q=S/180*(+f||0),c=[],t,x=a._.cacher(function(a,b,c){var d=a*F.cos(c)-b*F.sin(c);a=a*F.sin(c)+b*F.cos(c);return{x:d,y:a}});if(s)v=s[0],t=s[1],l=s[2],u=s[3];else{t=x(b,d,-q);b=t.x;d=t.y;t=x(k,p,-q);k=t.x;p=t.y;F.cos(S/180*f);F.sin(S/180*f);t=(b-k)/2;v=(d-p)/2;u=t*t/(e*e)+v*v/(h*h);1<u&&(u=F.sqrt(u),e*=u,h*=u);var u=e*e,w=h*h,u=(l==n?-1:1)*F.sqrt(Z((u*w-u*v*v-w*t*t)/
(u*v*v+w*t*t)));l=u*e*v/h+(b+k)/2;var u=u*-h*t/e+(d+p)/2,v=F.asin(((d-u)/h).toFixed(9));t=F.asin(((p-u)/h).toFixed(9));v=b<l?S-v:v;t=k<l?S-t:t;0>v&&(v=2*S+v);0>t&&(t=2*S+t);n&&v>t&&(v-=2*S);!n&&t>v&&(t-=2*S)}if(Z(t-v)>r){var c=t,w=k,G=p;t=v+r*(n&&t>v?1:-1);k=l+e*F.cos(t);p=u+h*F.sin(t);c=K(k,p,e,h,f,0,n,w,G,[t,c,l,u])}l=t-v;f=F.cos(v);r=F.sin(v);n=F.cos(t);t=F.sin(t);l=F.tan(l/4);e=4/3*e*l;l*=4/3*h;h=[b,d];b=[b+e*r,d-l*f];d=[k+e*t,p-l*n];k=[k,p];b[0]=2*h[0]-b[0];b[1]=2*h[1]-b[1];if(s)return[b,d,k].concat(c);
c=[b,d,k].concat(c).join().split(",");s=[];k=0;for(p=c.length;k<p;k++)s[k]=k%2?x(c[k-1],c[k],q).y:x(c[k],c[k+1],q).x;return s}function U(a,b,d,e,h,f,l,k){for(var n=[],p=[[],[] ],s,r,c,t,q=0;2>q;++q)0==q?(r=6*a-12*d+6*h,s=-3*a+9*d-9*h+3*l,c=3*d-3*a):(r=6*b-12*e+6*f,s=-3*b+9*e-9*f+3*k,c=3*e-3*b),1E-12>Z(s)?1E-12>Z(r)||(s=-c/r,0<s&&1>s&&n.push(s)):(t=r*r-4*c*s,c=F.sqrt(t),0>t||(t=(-r+c)/(2*s),0<t&&1>t&&n.push(t),s=(-r-c)/(2*s),0<s&&1>s&&n.push(s)));for(r=q=n.length;q--;)s=n[q],c=1-s,p[0][q]=c*c*c*a+3*
c*c*s*d+3*c*s*s*h+s*s*s*l,p[1][q]=c*c*c*b+3*c*c*s*e+3*c*s*s*f+s*s*s*k;p[0][r]=a;p[1][r]=b;p[0][r+1]=l;p[1][r+1]=k;p[0].length=p[1].length=r+2;return{min:{x:X.apply(0,p[0]),y:X.apply(0,p[1])},max:{x:W.apply(0,p[0]),y:W.apply(0,p[1])}}}function I(a,b){var e=!b&&A(a);if(!b&&e.curve)return d(e.curve);var f=G(a),l=b&&G(b),n={x:0,y:0,bx:0,by:0,X:0,Y:0,qx:null,qy:null},k={x:0,y:0,bx:0,by:0,X:0,Y:0,qx:null,qy:null},p=function(a,b,c){if(!a)return["C",b.x,b.y,b.x,b.y,b.x,b.y];a[0]in{T:1,Q:1}||(b.qx=b.qy=null);
switch(a[0]){case "M":b.X=a[1];b.Y=a[2];break;case "A":a=["C"].concat(K.apply(0,[b.x,b.y].concat(a.slice(1))));break;case "S":"C"==c||"S"==c?(c=2*b.x-b.bx,b=2*b.y-b.by):(c=b.x,b=b.y);a=["C",c,b].concat(a.slice(1));break;case "T":"Q"==c||"T"==c?(b.qx=2*b.x-b.qx,b.qy=2*b.y-b.qy):(b.qx=b.x,b.qy=b.y);a=["C"].concat(J(b.x,b.y,b.qx,b.qy,a[1],a[2]));break;case "Q":b.qx=a[1];b.qy=a[2];a=["C"].concat(J(b.x,b.y,a[1],a[2],a[3],a[4]));break;case "L":a=["C"].concat(h(b.x,b.y,a[1],a[2]));break;case "H":a=["C"].concat(h(b.x,
b.y,a[1],b.y));break;case "V":a=["C"].concat(h(b.x,b.y,b.x,a[1]));break;case "Z":a=["C"].concat(h(b.x,b.y,b.X,b.Y))}return a},s=function(a,b){if(7<a[b].length){a[b].shift();for(var c=a[b];c.length;)q[b]="A",l&&(u[b]="A"),a.splice(b++,0,["C"].concat(c.splice(0,6)));a.splice(b,1);v=W(f.length,l&&l.length||0)}},r=function(a,b,c,d,e){a&&b&&"M"==a[e][0]&&"M"!=b[e][0]&&(b.splice(e,0,["M",d.x,d.y]),c.bx=0,c.by=0,c.x=a[e][1],c.y=a[e][2],v=W(f.length,l&&l.length||0))},q=[],u=[],c="",t="",x=0,v=W(f.length,
l&&l.length||0);for(;x<v;x++){f[x]&&(c=f[x][0]);"C"!=c&&(q[x]=c,x&&(t=q[x-1]));f[x]=p(f[x],n,t);"A"!=q[x]&&"C"==c&&(q[x]="C");s(f,x);l&&(l[x]&&(c=l[x][0]),"C"!=c&&(u[x]=c,x&&(t=u[x-1])),l[x]=p(l[x],k,t),"A"!=u[x]&&"C"==c&&(u[x]="C"),s(l,x));r(f,l,n,k,x);r(l,f,k,n,x);var w=f[x],z=l&&l[x],y=w.length,U=l&&z.length;n.x=w[y-2];n.y=w[y-1];n.bx=$(w[y-4])||n.x;n.by=$(w[y-3])||n.y;k.bx=l&&($(z[U-4])||k.x);k.by=l&&($(z[U-3])||k.y);k.x=l&&z[U-2];k.y=l&&z[U-1]}l||(e.curve=d(f));return l?[f,l]:f}function P(a,
b){for(var d=[],e=0,h=a.length;h-2*!b>e;e+=2){var f=[{x:+a[e-2],y:+a[e-1]},{x:+a[e],y:+a[e+1]},{x:+a[e+2],y:+a[e+3]},{x:+a[e+4],y:+a[e+5]}];b?e?h-4==e?f[3]={x:+a[0],y:+a[1]}:h-2==e&&(f[2]={x:+a[0],y:+a[1]},f[3]={x:+a[2],y:+a[3]}):f[0]={x:+a[h-2],y:+a[h-1]}:h-4==e?f[3]=f[2]:e||(f[0]={x:+a[e],y:+a[e+1]});d.push(["C",(-f[0].x+6*f[1].x+f[2].x)/6,(-f[0].y+6*f[1].y+f[2].y)/6,(f[1].x+6*f[2].x-f[3].x)/6,(f[1].y+6*f[2].y-f[3].y)/6,f[2].x,f[2].y])}return d}y=k.prototype;var Q=a.is,C=a._.clone,L="hasOwnProperty",
N=/,?([a-z]),?/gi,$=parseFloat,F=Math,S=F.PI,X=F.min,W=F.max,ma=F.pow,Z=F.abs;M=n(1);var na=n(),ba=n(0,1),V=a._unit2px;a.path=A;a.path.getTotalLength=M;a.path.getPointAtLength=na;a.path.getSubpath=function(a,b,d){if(1E-6>this.getTotalLength(a)-d)return ba(a,b).end;a=ba(a,d,1);return b?ba(a,b).end:a};y.getTotalLength=function(){if(this.node.getTotalLength)return this.node.getTotalLength()};y.getPointAtLength=function(a){return na(this.attr("d"),a)};y.getSubpath=function(b,d){return a.path.getSubpath(this.attr("d"),
b,d)};a._.box=w;a.path.findDotsAtSegment=u;a.path.bezierBBox=p;a.path.isPointInsideBBox=b;a.path.isBBoxIntersect=q;a.path.intersection=function(a,b){return l(a,b)};a.path.intersectionNumber=function(a,b){return l(a,b,1)};a.path.isPointInside=function(a,d,e){var h=r(a);return b(h,d,e)&&1==l(a,[["M",d,e],["H",h.x2+10] ],1)%2};a.path.getBBox=r;a.path.get={path:function(a){return a.attr("path")},circle:function(a){a=V(a);return x(a.cx,a.cy,a.r)},ellipse:function(a){a=V(a);return x(a.cx||0,a.cy||0,a.rx,
a.ry)},rect:function(a){a=V(a);return s(a.x||0,a.y||0,a.width,a.height,a.rx,a.ry)},image:function(a){a=V(a);return s(a.x||0,a.y||0,a.width,a.height)},line:function(a){return"M"+[a.attr("x1")||0,a.attr("y1")||0,a.attr("x2"),a.attr("y2")]},polyline:function(a){return"M"+a.attr("points")},polygon:function(a){return"M"+a.attr("points")+"z"},deflt:function(a){a=a.node.getBBox();return s(a.x,a.y,a.width,a.height)}};a.path.toRelative=function(b){var e=A(b),h=String.prototype.toLowerCase;if(e.rel)return d(e.rel);
a.is(b,"array")&&a.is(b&&b[0],"array")||(b=a.parsePathString(b));var f=[],l=0,n=0,k=0,p=0,s=0;"M"==b[0][0]&&(l=b[0][1],n=b[0][2],k=l,p=n,s++,f.push(["M",l,n]));for(var r=b.length;s<r;s++){var q=f[s]=[],x=b[s];if(x[0]!=h.call(x[0]))switch(q[0]=h.call(x[0]),q[0]){case "a":q[1]=x[1];q[2]=x[2];q[3]=x[3];q[4]=x[4];q[5]=x[5];q[6]=+(x[6]-l).toFixed(3);q[7]=+(x[7]-n).toFixed(3);break;case "v":q[1]=+(x[1]-n).toFixed(3);break;case "m":k=x[1],p=x[2];default:for(var c=1,t=x.length;c<t;c++)q[c]=+(x[c]-(c%2?l:
n)).toFixed(3)}else for(f[s]=[],"m"==x[0]&&(k=x[1]+l,p=x[2]+n),q=0,c=x.length;q<c;q++)f[s][q]=x[q];x=f[s].length;switch(f[s][0]){case "z":l=k;n=p;break;case "h":l+=+f[s][x-1];break;case "v":n+=+f[s][x-1];break;default:l+=+f[s][x-2],n+=+f[s][x-1]}}f.toString=z;e.rel=d(f);return f};a.path.toAbsolute=G;a.path.toCubic=I;a.path.map=function(a,b){if(!b)return a;var d,e,h,f,l,n,k;a=I(a);h=0;for(l=a.length;h<l;h++)for(k=a[h],f=1,n=k.length;f<n;f+=2)d=b.x(k[f],k[f+1]),e=b.y(k[f],k[f+1]),k[f]=d,k[f+1]=e;return a};
a.path.toString=z;a.path.clone=d});C.plugin(function(a,v,y,C){var A=Math.max,w=Math.min,z=function(a){this.items=[];this.bindings={};this.length=0;this.type="set";if(a)for(var f=0,n=a.length;f<n;f++)a[f]&&(this[this.items.length]=this.items[this.items.length]=a[f],this.length++)};v=z.prototype;v.push=function(){for(var a,f,n=0,k=arguments.length;n<k;n++)if(a=arguments[n])f=this.items.length,this[f]=this.items[f]=a,this.length++;return this};v.pop=function(){this.length&&delete this[this.length--];
return this.items.pop()};v.forEach=function(a,f){for(var n=0,k=this.items.length;n<k&&!1!==a.call(f,this.items[n],n);n++);return this};v.animate=function(d,f,n,u){"function"!=typeof n||n.length||(u=n,n=L.linear);d instanceof a._.Animation&&(u=d.callback,n=d.easing,f=n.dur,d=d.attr);var p=arguments;if(a.is(d,"array")&&a.is(p[p.length-1],"array"))var b=!0;var q,e=function(){q?this.b=q:q=this.b},l=0,r=u&&function(){l++==this.length&&u.call(this)};return this.forEach(function(a,l){k.once("snap.animcreated."+
a.id,e);b?p[l]&&a.animate.apply(a,p[l]):a.animate(d,f,n,r)})};v.remove=function(){for(;this.length;)this.pop().remove();return this};v.bind=function(a,f,k){var u={};if("function"==typeof f)this.bindings[a]=f;else{var p=k||a;this.bindings[a]=function(a){u[p]=a;f.attr(u)}}return this};v.attr=function(a){var f={},k;for(k in a)if(this.bindings[k])this.bindings[k](a[k]);else f[k]=a[k];a=0;for(k=this.items.length;a<k;a++)this.items[a].attr(f);return this};v.clear=function(){for(;this.length;)this.pop()};
v.splice=function(a,f,k){a=0>a?A(this.length+a,0):a;f=A(0,w(this.length-a,f));var u=[],p=[],b=[],q;for(q=2;q<arguments.length;q++)b.push(arguments[q]);for(q=0;q<f;q++)p.push(this[a+q]);for(;q<this.length-a;q++)u.push(this[a+q]);var e=b.length;for(q=0;q<e+u.length;q++)this.items[a+q]=this[a+q]=q<e?b[q]:u[q-e];for(q=this.items.length=this.length-=f-e;this[q];)delete this[q++];return new z(p)};v.exclude=function(a){for(var f=0,k=this.length;f<k;f++)if(this[f]==a)return this.splice(f,1),!0;return!1};
v.insertAfter=function(a){for(var f=this.items.length;f--;)this.items[f].insertAfter(a);return this};v.getBBox=function(){for(var a=[],f=[],k=[],u=[],p=this.items.length;p--;)if(!this.items[p].removed){var b=this.items[p].getBBox();a.push(b.x);f.push(b.y);k.push(b.x+b.width);u.push(b.y+b.height)}a=w.apply(0,a);f=w.apply(0,f);k=A.apply(0,k);u=A.apply(0,u);return{x:a,y:f,x2:k,y2:u,width:k-a,height:u-f,cx:a+(k-a)/2,cy:f+(u-f)/2}};v.clone=function(a){a=new z;for(var f=0,k=this.items.length;f<k;f++)a.push(this.items[f].clone());
return a};v.toString=function(){return"Snap\u2018s set"};v.type="set";a.set=function(){var a=new z;arguments.length&&a.push.apply(a,Array.prototype.slice.call(arguments,0));return a}});C.plugin(function(a,v,y,C){function A(a){var b=a[0];switch(b.toLowerCase()){case "t":return[b,0,0];case "m":return[b,1,0,0,1,0,0];case "r":return 4==a.length?[b,0,a[2],a[3] ]:[b,0];case "s":return 5==a.length?[b,1,1,a[3],a[4] ]:3==a.length?[b,1,1]:[b,1]}}function w(b,d,f){d=q(d).replace(/\.{3}|\u2026/g,b);b=a.parseTransformString(b)||
[];d=a.parseTransformString(d)||[];for(var k=Math.max(b.length,d.length),p=[],v=[],h=0,w,z,y,I;h<k;h++){y=b[h]||A(d[h]);I=d[h]||A(y);if(y[0]!=I[0]||"r"==y[0].toLowerCase()&&(y[2]!=I[2]||y[3]!=I[3])||"s"==y[0].toLowerCase()&&(y[3]!=I[3]||y[4]!=I[4])){b=a._.transform2matrix(b,f());d=a._.transform2matrix(d,f());p=[["m",b.a,b.b,b.c,b.d,b.e,b.f] ];v=[["m",d.a,d.b,d.c,d.d,d.e,d.f] ];break}p[h]=[];v[h]=[];w=0;for(z=Math.max(y.length,I.length);w<z;w++)w in y&&(p[h][w]=y[w]),w in I&&(v[h][w]=I[w])}return{from:u(p),
to:u(v),f:n(p)}}function z(a){return a}function d(a){return function(b){return+b.toFixed(3)+a}}function f(b){return a.rgb(b[0],b[1],b[2])}function n(a){var b=0,d,f,k,n,h,p,q=[];d=0;for(f=a.length;d<f;d++){h="[";p=['"'+a[d][0]+'"'];k=1;for(n=a[d].length;k<n;k++)p[k]="val["+b++ +"]";h+=p+"]";q[d]=h}return Function("val","return Snap.path.toString.call(["+q+"])")}function u(a){for(var b=[],d=0,f=a.length;d<f;d++)for(var k=1,n=a[d].length;k<n;k++)b.push(a[d][k]);return b}var p={},b=/[a-z]+$/i,q=String;
p.stroke=p.fill="colour";v.prototype.equal=function(a,b){return k("snap.util.equal",this,a,b).firstDefined()};k.on("snap.util.equal",function(e,k){var r,s;r=q(this.attr(e)||"");var x=this;if(r==+r&&k==+k)return{from:+r,to:+k,f:z};if("colour"==p[e])return r=a.color(r),s=a.color(k),{from:[r.r,r.g,r.b,r.opacity],to:[s.r,s.g,s.b,s.opacity],f:f};if("transform"==e||"gradientTransform"==e||"patternTransform"==e)return k instanceof a.Matrix&&(k=k.toTransformString()),a._.rgTransform.test(k)||(k=a._.svgTransform2string(k)),
w(r,k,function(){return x.getBBox(1)});if("d"==e||"path"==e)return r=a.path.toCubic(r,k),{from:u(r[0]),to:u(r[1]),f:n(r[0])};if("points"==e)return r=q(r).split(a._.separator),s=q(k).split(a._.separator),{from:r,to:s,f:function(a){return a}};aUnit=r.match(b);s=q(k).match(b);return aUnit&&aUnit==s?{from:parseFloat(r),to:parseFloat(k),f:d(aUnit)}:{from:this.asPX(e),to:this.asPX(e,k),f:z}})});C.plugin(function(a,v,y,C){var A=v.prototype,w="createTouch"in C.doc;v="click dblclick mousedown mousemove mouseout mouseover mouseup touchstart touchmove touchend touchcancel".split(" ");
var z={mousedown:"touchstart",mousemove:"touchmove",mouseup:"touchend"},d=function(a,b){var d="y"==a?"scrollTop":"scrollLeft",e=b&&b.node?b.node.ownerDocument:C.doc;return e[d in e.documentElement?"documentElement":"body"][d]},f=function(){this.returnValue=!1},n=function(){return this.originalEvent.preventDefault()},u=function(){this.cancelBubble=!0},p=function(){return this.originalEvent.stopPropagation()},b=function(){if(C.doc.addEventListener)return function(a,b,e,f){var k=w&&z[b]?z[b]:b,l=function(k){var l=
d("y",f),q=d("x",f);if(w&&z.hasOwnProperty(b))for(var r=0,u=k.targetTouches&&k.targetTouches.length;r<u;r++)if(k.targetTouches[r].target==a||a.contains(k.targetTouches[r].target)){u=k;k=k.targetTouches[r];k.originalEvent=u;k.preventDefault=n;k.stopPropagation=p;break}return e.call(f,k,k.clientX+q,k.clientY+l)};b!==k&&a.addEventListener(b,l,!1);a.addEventListener(k,l,!1);return function(){b!==k&&a.removeEventListener(b,l,!1);a.removeEventListener(k,l,!1);return!0}};if(C.doc.attachEvent)return function(a,
b,e,h){var k=function(a){a=a||h.node.ownerDocument.window.event;var b=d("y",h),k=d("x",h),k=a.clientX+k,b=a.clientY+b;a.preventDefault=a.preventDefault||f;a.stopPropagation=a.stopPropagation||u;return e.call(h,a,k,b)};a.attachEvent("on"+b,k);return function(){a.detachEvent("on"+b,k);return!0}}}(),q=[],e=function(a){for(var b=a.clientX,e=a.clientY,f=d("y"),l=d("x"),n,p=q.length;p--;){n=q[p];if(w)for(var r=a.touches&&a.touches.length,u;r--;){if(u=a.touches[r],u.identifier==n.el._drag.id||n.el.node.contains(u.target)){b=
u.clientX;e=u.clientY;(a.originalEvent?a.originalEvent:a).preventDefault();break}}else a.preventDefault();b+=l;e+=f;k("snap.drag.move."+n.el.id,n.move_scope||n.el,b-n.el._drag.x,e-n.el._drag.y,b,e,a)}},l=function(b){a.unmousemove(e).unmouseup(l);for(var d=q.length,f;d--;)f=q[d],f.el._drag={},k("snap.drag.end."+f.el.id,f.end_scope||f.start_scope||f.move_scope||f.el,b);q=[]};for(y=v.length;y--;)(function(d){a[d]=A[d]=function(e,f){a.is(e,"function")&&(this.events=this.events||[],this.events.push({name:d,
f:e,unbind:b(this.node||document,d,e,f||this)}));return this};a["un"+d]=A["un"+d]=function(a){for(var b=this.events||[],e=b.length;e--;)if(b[e].name==d&&(b[e].f==a||!a)){b[e].unbind();b.splice(e,1);!b.length&&delete this.events;break}return this}})(v[y]);A.hover=function(a,b,d,e){return this.mouseover(a,d).mouseout(b,e||d)};A.unhover=function(a,b){return this.unmouseover(a).unmouseout(b)};var r=[];A.drag=function(b,d,f,h,n,p){function u(r,v,w){(r.originalEvent||r).preventDefault();this._drag.x=v;
this._drag.y=w;this._drag.id=r.identifier;!q.length&&a.mousemove(e).mouseup(l);q.push({el:this,move_scope:h,start_scope:n,end_scope:p});d&&k.on("snap.drag.start."+this.id,d);b&&k.on("snap.drag.move."+this.id,b);f&&k.on("snap.drag.end."+this.id,f);k("snap.drag.start."+this.id,n||h||this,v,w,r)}if(!arguments.length){var v;return this.drag(function(a,b){this.attr({transform:v+(v?"T":"t")+[a,b]})},function(){v=this.transform().local})}this._drag={};r.push({el:this,start:u});this.mousedown(u);return this};
A.undrag=function(){for(var b=r.length;b--;)r[b].el==this&&(this.unmousedown(r[b].start),r.splice(b,1),k.unbind("snap.drag.*."+this.id));!r.length&&a.unmousemove(e).unmouseup(l);return this}});C.plugin(function(a,v,y,C){y=y.prototype;var A=/^\s*url\((.+)\)/,w=String,z=a._.$;a.filter={};y.filter=function(d){var f=this;"svg"!=f.type&&(f=f.paper);d=a.parse(w(d));var k=a._.id(),u=z("filter");z(u,{id:k,filterUnits:"userSpaceOnUse"});u.appendChild(d.node);f.defs.appendChild(u);return new v(u)};k.on("snap.util.getattr.filter",
function(){k.stop();var d=z(this.node,"filter");if(d)return(d=w(d).match(A))&&a.select(d[1])});k.on("snap.util.attr.filter",function(d){if(d instanceof v&&"filter"==d.type){k.stop();var f=d.node.id;f||(z(d.node,{id:d.id}),f=d.id);z(this.node,{filter:a.url(f)})}d&&"none"!=d||(k.stop(),this.node.removeAttribute("filter"))});a.filter.blur=function(d,f){null==d&&(d=2);return a.format('<feGaussianBlur stdDeviation="{def}"/>',{def:null==f?d:[d,f]})};a.filter.blur.toString=function(){return this()};a.filter.shadow=
function(d,f,k,u,p){"string"==typeof k&&(p=u=k,k=4);"string"!=typeof u&&(p=u,u="#000");null==k&&(k=4);null==p&&(p=1);null==d&&(d=0,f=2);null==f&&(f=d);u=a.color(u||"#000");return a.format('<feGaussianBlur in="SourceAlpha" stdDeviation="{blur}"/><feOffset dx="{dx}" dy="{dy}" result="offsetblur"/><feFlood flood-color="{color}"/><feComposite in2="offsetblur" operator="in"/><feComponentTransfer><feFuncA type="linear" slope="{opacity}"/></feComponentTransfer><feMerge><feMergeNode/><feMergeNode in="SourceGraphic"/></feMerge>',
{color:u,dx:d,dy:f,blur:k,opacity:p})};a.filter.shadow.toString=function(){return this()};a.filter.grayscale=function(d){null==d&&(d=1);return a.format('<feColorMatrix type="matrix" values="{a} {b} {c} 0 0 {d} {e} {f} 0 0 {g} {b} {h} 0 0 0 0 0 1 0"/>',{a:0.2126+0.7874*(1-d),b:0.7152-0.7152*(1-d),c:0.0722-0.0722*(1-d),d:0.2126-0.2126*(1-d),e:0.7152+0.2848*(1-d),f:0.0722-0.0722*(1-d),g:0.2126-0.2126*(1-d),h:0.0722+0.9278*(1-d)})};a.filter.grayscale.toString=function(){return this()};a.filter.sepia=
function(d){null==d&&(d=1);return a.format('<feColorMatrix type="matrix" values="{a} {b} {c} 0 0 {d} {e} {f} 0 0 {g} {h} {i} 0 0 0 0 0 1 0"/>',{a:0.393+0.607*(1-d),b:0.769-0.769*(1-d),c:0.189-0.189*(1-d),d:0.349-0.349*(1-d),e:0.686+0.314*(1-d),f:0.168-0.168*(1-d),g:0.272-0.272*(1-d),h:0.534-0.534*(1-d),i:0.131+0.869*(1-d)})};a.filter.sepia.toString=function(){return this()};a.filter.saturate=function(d){null==d&&(d=1);return a.format('<feColorMatrix type="saturate" values="{amount}"/>',{amount:1-
d})};a.filter.saturate.toString=function(){return this()};a.filter.hueRotate=function(d){return a.format('<feColorMatrix type="hueRotate" values="{angle}"/>',{angle:d||0})};a.filter.hueRotate.toString=function(){return this()};a.filter.invert=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="table" tableValues="{amount} {amount2}"/><feFuncG type="table" tableValues="{amount} {amount2}"/><feFuncB type="table" tableValues="{amount} {amount2}"/></feComponentTransfer>',{amount:d,
amount2:1-d})};a.filter.invert.toString=function(){return this()};a.filter.brightness=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="linear" slope="{amount}"/><feFuncG type="linear" slope="{amount}"/><feFuncB type="linear" slope="{amount}"/></feComponentTransfer>',{amount:d})};a.filter.brightness.toString=function(){return this()};a.filter.contrast=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="linear" slope="{amount}" intercept="{amount2}"/><feFuncG type="linear" slope="{amount}" intercept="{amount2}"/><feFuncB type="linear" slope="{amount}" intercept="{amount2}"/></feComponentTransfer>',
{amount:d,amount2:0.5-d/2})};a.filter.contrast.toString=function(){return this()}});return C});

]]> </script>
</svg>


<pre class='language-julia'><code class='language-julia'>mat(put(3, 1=&gt;H))</code></pre>
<pre id='var-hash169404' class='code-output documenter-example-output'>8×8 SparseMatrixCSC{ComplexF64, Int64} with 16 stored entries:
 0.707107+0.0im   0.707107+0.0im           ⋅      …           ⋅                ⋅    
 0.707107+0.0im  -0.707107+0.0im           ⋅                  ⋅                ⋅    
          ⋅                ⋅      0.707107+0.0im              ⋅                ⋅    
          ⋅                ⋅      0.707107+0.0im              ⋅                ⋅    
          ⋅                ⋅               ⋅                  ⋅                ⋅    
          ⋅                ⋅               ⋅      …           ⋅                ⋅    
          ⋅                ⋅               ⋅         0.707107+0.0im   0.707107+0.0im
          ⋅                ⋅               ⋅         0.707107+0.0im  -0.707107+0.0im</pre>


<div class="markdown"><p class="tex">$$&#123;\rm shift&#125;&#40;\theta&#41; &#61; \begin&#123;pmatrix&#125;
1 &amp; 0\\
0 &amp; e^&#123;i\theta&#125;
\end&#123;pmatrix&#125;$$</p>
</div>

<pre class='language-julia'><code class='language-julia'>mat(shift(π/2))</code></pre>
<pre id='var-hash248302' class='code-output documenter-example-output'>2×2 Diagonal{ComplexF64, Vector{ComplexF64}}:
 1.0+0.0im              ⋅    
     ⋅      6.12323e-17+1.0im</pre>


<div class="markdown"><p>Control block is a composite block, its matrix representation is</p>
</div>


<div class="markdown"><p class="tex">$$\begin&#123;align&#125;
\texttt&#123;control&#40;n, i, j&#61;&gt;G&#41;&#125; &#61; &amp;I^&#123;\otimes n-i&#125; P_0 \otimes I^&#123;\otimes i-j-1&#125; \otimes I\otimes I^&#123;\otimes j-1&#125;
&#43;\\
&amp; I^&#123;\otimes n-i&#125; P_1 \otimes I^&#123;\otimes i-j-1&#125; \otimes G\otimes I^&#123;\otimes j-1&#125;
\end&#123;align&#125;$$</p>
</div>


<div class="markdown"><p><span class="tex">$P_0&#61;|0\rangle\langle 0|$</span></p>
</div>

<pre class='language-julia'><code class='language-julia'>mat(ConstGate.P0)</code></pre>
<pre id='var-hash142206' class='code-output documenter-example-output'>2×2 Diagonal{ComplexF64, Vector{ComplexF64}}:
 1.0+0.0im      ⋅    
     ⋅      0.0+0.0im</pre>


<div class="markdown"><p><span class="tex">$P_1&#61;|1\rangle\langle 1|$</span></p>
</div>

<pre class='language-julia'><code class='language-julia'>mat(ConstGate.P1)</code></pre>
<pre id='var-hash700734' class='code-output documenter-example-output'>2×2 Diagonal{ComplexF64, Vector{ComplexF64}}:
 0.0+0.0im      ⋅    
     ⋅      1.0+0.0im</pre>

<pre class='language-julia'><code class='language-julia'>vizcircuit(control(3, 3, 1=&gt;shift(π/4)))</code></pre>
<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:xlink="http://www.w3.org/1999/xlink"
     version="1.2"
     width="42.5mm" height="56.25mm" viewBox="0 0 42.5 56.25"
     stroke="none"
     fill="#000000"
     stroke-width="0.3"
     font-size="3.88"

     id="img-da86ccac">
<defs>
  <marker id="arrow" markerWidth="15" markerHeight="7" refX="5" refY="3.5" orient="auto" markerUnits="strokeWidth">
    <path d="M0,0 L15,3.5 L0,7 z" stroke="context-stroke" fill="context-stroke"/>
  </marker>
</defs>
<g stroke-width="0" fill="#000000" id="img-da86ccac-1">
  <g transform="translate(23.12,43.75)">
    <circle cx="0" cy="0" r="0.83" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-da86ccac-2">
  <g transform="translate(23.12,6.25)">
    <path d="M-7.2,-5 L7.2,-5 7.2,5 -7.2,5  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-da86ccac-3">
  <g transform="translate(33.75,43.75)">
    <path fill="none" d="M0,0 L0,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-da86ccac-4">
  <g transform="translate(33.75,25)">
    <path fill="none" d="M0,0 L0,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-da86ccac-5">
  <g transform="translate(33.75,6.25)">
    <path fill="none" d="M0,0 L0,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-da86ccac-6">
  <g transform="translate(23.12,25)">
    <path fill="none" d="M-10.62,0 L10.62,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-da86ccac-7">
  <g transform="translate(28.85,43.75)">
    <path fill="none" d="M-4.9,0 L4.9,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-da86ccac-8">
  <g transform="translate(17.4,43.75)">
    <path fill="none" d="M-4.9,0 L4.9,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-da86ccac-9">
  <g transform="translate(32.04,6.25)">
    <path fill="none" d="M-1.71,0 L1.71,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-da86ccac-10">
  <g transform="translate(14.21,6.25)">
    <path fill="none" d="M-1.71,0 L1.71,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-da86ccac-11">
  <g transform="translate(23.12,27.08)">
    <path fill="none" d="M0,15.83 L0,-15.83 " class="primitive"/>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="5.64" id="img-da86ccac-12">
  <g transform="translate(23.12,43.75)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em"></text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="3.53" id="img-da86ccac-13">
  <g transform="translate(23.12,6.25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">ϕ(π/4)</text>
    </g>
  </g>
</g>
<script> <![CDATA[
(function(N){var k=/[\.\/]/,L=/\s*,\s*/,C=function(a,d){return a-d},a,v,y={n:{}},M=function(){for(var a=0,d=this.length;a<d;a++)if("undefined"!=typeof this[a])return this[a]},A=function(){for(var a=this.length;--a;)if("undefined"!=typeof this[a])return this[a]},w=function(k,d){k=String(k);var f=v,n=Array.prototype.slice.call(arguments,2),u=w.listeners(k),p=0,b,q=[],e={},l=[],r=a;l.firstDefined=M;l.lastDefined=A;a=k;for(var s=v=0,x=u.length;s<x;s++)"zIndex"in u[s]&&(q.push(u[s].zIndex),0>u[s].zIndex&&
(e[u[s].zIndex]=u[s]));for(q.sort(C);0>q[p];)if(b=e[q[p++] ],l.push(b.apply(d,n)),v)return v=f,l;for(s=0;s<x;s++)if(b=u[s],"zIndex"in b)if(b.zIndex==q[p]){l.push(b.apply(d,n));if(v)break;do if(p++,(b=e[q[p] ])&&l.push(b.apply(d,n)),v)break;while(b)}else e[b.zIndex]=b;else if(l.push(b.apply(d,n)),v)break;v=f;a=r;return l};w._events=y;w.listeners=function(a){a=a.split(k);var d=y,f,n,u,p,b,q,e,l=[d],r=[];u=0;for(p=a.length;u<p;u++){e=[];b=0;for(q=l.length;b<q;b++)for(d=l[b].n,f=[d[a[u] ],d["*"] ],n=2;n--;)if(d=
f[n])e.push(d),r=r.concat(d.f||[]);l=e}return r};w.on=function(a,d){a=String(a);if("function"!=typeof d)return function(){};for(var f=a.split(L),n=0,u=f.length;n<u;n++)(function(a){a=a.split(k);for(var b=y,f,e=0,l=a.length;e<l;e++)b=b.n,b=b.hasOwnProperty(a[e])&&b[a[e] ]||(b[a[e] ]={n:{}});b.f=b.f||[];e=0;for(l=b.f.length;e<l;e++)if(b.f[e]==d){f=!0;break}!f&&b.f.push(d)})(f[n]);return function(a){+a==+a&&(d.zIndex=+a)}};w.f=function(a){var d=[].slice.call(arguments,1);return function(){w.apply(null,
[a,null].concat(d).concat([].slice.call(arguments,0)))}};w.stop=function(){v=1};w.nt=function(k){return k?(new RegExp("(?:\\.|\\/|^)"+k+"(?:\\.|\\/|$)")).test(a):a};w.nts=function(){return a.split(k)};w.off=w.unbind=function(a,d){if(a){var f=a.split(L);if(1<f.length)for(var n=0,u=f.length;n<u;n++)w.off(f[n],d);else{for(var f=a.split(k),p,b,q,e,l=[y],n=0,u=f.length;n<u;n++)for(e=0;e<l.length;e+=q.length-2){q=[e,1];p=l[e].n;if("*"!=f[n])p[f[n] ]&&q.push(p[f[n] ]);else for(b in p)p.hasOwnProperty(b)&&
q.push(p[b]);l.splice.apply(l,q)}n=0;for(u=l.length;n<u;n++)for(p=l[n];p.n;){if(d){if(p.f){e=0;for(f=p.f.length;e<f;e++)if(p.f[e]==d){p.f.splice(e,1);break}!p.f.length&&delete p.f}for(b in p.n)if(p.n.hasOwnProperty(b)&&p.n[b].f){q=p.n[b].f;e=0;for(f=q.length;e<f;e++)if(q[e]==d){q.splice(e,1);break}!q.length&&delete p.n[b].f}}else for(b in delete p.f,p.n)p.n.hasOwnProperty(b)&&p.n[b].f&&delete p.n[b].f;p=p.n}}}else w._events=y={n:{}}};w.once=function(a,d){var f=function(){w.unbind(a,f);return d.apply(this,
arguments)};return w.on(a,f)};w.version="0.4.2";w.toString=function(){return"You are running Eve 0.4.2"};"undefined"!=typeof module&&module.exports?module.exports=w:"function"===typeof define&&define.amd?define("eve",[],function(){return w}):N.eve=w})(this);
(function(N,k){"function"===typeof define&&define.amd?define("Snap.svg",["eve"],function(L){return k(N,L)}):k(N,N.eve)})(this,function(N,k){var L=function(a){var k={},y=N.requestAnimationFrame||N.webkitRequestAnimationFrame||N.mozRequestAnimationFrame||N.oRequestAnimationFrame||N.msRequestAnimationFrame||function(a){setTimeout(a,16)},M=Array.isArray||function(a){return a instanceof Array||"[object Array]"==Object.prototype.toString.call(a)},A=0,w="M"+(+new Date).toString(36),z=function(a){if(null==
a)return this.s;var b=this.s-a;this.b+=this.dur*b;this.B+=this.dur*b;this.s=a},d=function(a){if(null==a)return this.spd;this.spd=a},f=function(a){if(null==a)return this.dur;this.s=this.s*a/this.dur;this.dur=a},n=function(){delete k[this.id];this.update();a("mina.stop."+this.id,this)},u=function(){this.pdif||(delete k[this.id],this.update(),this.pdif=this.get()-this.b)},p=function(){this.pdif&&(this.b=this.get()-this.pdif,delete this.pdif,k[this.id]=this)},b=function(){var a;if(M(this.start)){a=[];
for(var b=0,e=this.start.length;b<e;b++)a[b]=+this.start[b]+(this.end[b]-this.start[b])*this.easing(this.s)}else a=+this.start+(this.end-this.start)*this.easing(this.s);this.set(a)},q=function(){var l=0,b;for(b in k)if(k.hasOwnProperty(b)){var e=k[b],f=e.get();l++;e.s=(f-e.b)/(e.dur/e.spd);1<=e.s&&(delete k[b],e.s=1,l--,function(b){setTimeout(function(){a("mina.finish."+b.id,b)})}(e));e.update()}l&&y(q)},e=function(a,r,s,x,G,h,J){a={id:w+(A++).toString(36),start:a,end:r,b:s,s:0,dur:x-s,spd:1,get:G,
set:h,easing:J||e.linear,status:z,speed:d,duration:f,stop:n,pause:u,resume:p,update:b};k[a.id]=a;r=0;for(var K in k)if(k.hasOwnProperty(K)&&(r++,2==r))break;1==r&&y(q);return a};e.time=Date.now||function(){return+new Date};e.getById=function(a){return k[a]||null};e.linear=function(a){return a};e.easeout=function(a){return Math.pow(a,1.7)};e.easein=function(a){return Math.pow(a,0.48)};e.easeinout=function(a){if(1==a)return 1;if(0==a)return 0;var b=0.48-a/1.04,e=Math.sqrt(0.1734+b*b);a=e-b;a=Math.pow(Math.abs(a),
1/3)*(0>a?-1:1);b=-e-b;b=Math.pow(Math.abs(b),1/3)*(0>b?-1:1);a=a+b+0.5;return 3*(1-a)*a*a+a*a*a};e.backin=function(a){return 1==a?1:a*a*(2.70158*a-1.70158)};e.backout=function(a){if(0==a)return 0;a-=1;return a*a*(2.70158*a+1.70158)+1};e.elastic=function(a){return a==!!a?a:Math.pow(2,-10*a)*Math.sin(2*(a-0.075)*Math.PI/0.3)+1};e.bounce=function(a){a<1/2.75?a*=7.5625*a:a<2/2.75?(a-=1.5/2.75,a=7.5625*a*a+0.75):a<2.5/2.75?(a-=2.25/2.75,a=7.5625*a*a+0.9375):(a-=2.625/2.75,a=7.5625*a*a+0.984375);return a};
return N.mina=e}("undefined"==typeof k?function(){}:k),C=function(){function a(c,t){if(c){if(c.tagName)return x(c);if(y(c,"array")&&a.set)return a.set.apply(a,c);if(c instanceof e)return c;if(null==t)return c=G.doc.querySelector(c),x(c)}return new s(null==c?"100%":c,null==t?"100%":t)}function v(c,a){if(a){"#text"==c&&(c=G.doc.createTextNode(a.text||""));"string"==typeof c&&(c=v(c));if("string"==typeof a)return"xlink:"==a.substring(0,6)?c.getAttributeNS(m,a.substring(6)):"xml:"==a.substring(0,4)?c.getAttributeNS(la,
a.substring(4)):c.getAttribute(a);for(var da in a)if(a[h](da)){var b=J(a[da]);b?"xlink:"==da.substring(0,6)?c.setAttributeNS(m,da.substring(6),b):"xml:"==da.substring(0,4)?c.setAttributeNS(la,da.substring(4),b):c.setAttribute(da,b):c.removeAttribute(da)}}else c=G.doc.createElementNS(la,c);return c}function y(c,a){a=J.prototype.toLowerCase.call(a);return"finite"==a?isFinite(c):"array"==a&&(c instanceof Array||Array.isArray&&Array.isArray(c))?!0:"null"==a&&null===c||a==typeof c&&null!==c||"object"==
a&&c===Object(c)||$.call(c).slice(8,-1).toLowerCase()==a}function M(c){if("function"==typeof c||Object(c)!==c)return c;var a=new c.constructor,b;for(b in c)c[h](b)&&(a[b]=M(c[b]));return a}function A(c,a,b){function m(){var e=Array.prototype.slice.call(arguments,0),f=e.join("\u2400"),d=m.cache=m.cache||{},l=m.count=m.count||[];if(d[h](f)){a:for(var e=l,l=f,B=0,H=e.length;B<H;B++)if(e[B]===l){e.push(e.splice(B,1)[0]);break a}return b?b(d[f]):d[f]}1E3<=l.length&&delete d[l.shift()];l.push(f);d[f]=c.apply(a,
e);return b?b(d[f]):d[f]}return m}function w(c,a,b,m,e,f){return null==e?(c-=b,a-=m,c||a?(180*I.atan2(-a,-c)/C+540)%360:0):w(c,a,e,f)-w(b,m,e,f)}function z(c){return c%360*C/180}function d(c){var a=[];c=c.replace(/(?:^|\s)(\w+)\(([^)]+)\)/g,function(c,b,m){m=m.split(/\s*,\s*|\s+/);"rotate"==b&&1==m.length&&m.push(0,0);"scale"==b&&(2<m.length?m=m.slice(0,2):2==m.length&&m.push(0,0),1==m.length&&m.push(m[0],0,0));"skewX"==b?a.push(["m",1,0,I.tan(z(m[0])),1,0,0]):"skewY"==b?a.push(["m",1,I.tan(z(m[0])),
0,1,0,0]):a.push([b.charAt(0)].concat(m));return c});return a}function f(c,t){var b=O(c),m=new a.Matrix;if(b)for(var e=0,f=b.length;e<f;e++){var h=b[e],d=h.length,B=J(h[0]).toLowerCase(),H=h[0]!=B,l=H?m.invert():0,E;"t"==B&&2==d?m.translate(h[1],0):"t"==B&&3==d?H?(d=l.x(0,0),B=l.y(0,0),H=l.x(h[1],h[2]),l=l.y(h[1],h[2]),m.translate(H-d,l-B)):m.translate(h[1],h[2]):"r"==B?2==d?(E=E||t,m.rotate(h[1],E.x+E.width/2,E.y+E.height/2)):4==d&&(H?(H=l.x(h[2],h[3]),l=l.y(h[2],h[3]),m.rotate(h[1],H,l)):m.rotate(h[1],
h[2],h[3])):"s"==B?2==d||3==d?(E=E||t,m.scale(h[1],h[d-1],E.x+E.width/2,E.y+E.height/2)):4==d?H?(H=l.x(h[2],h[3]),l=l.y(h[2],h[3]),m.scale(h[1],h[1],H,l)):m.scale(h[1],h[1],h[2],h[3]):5==d&&(H?(H=l.x(h[3],h[4]),l=l.y(h[3],h[4]),m.scale(h[1],h[2],H,l)):m.scale(h[1],h[2],h[3],h[4])):"m"==B&&7==d&&m.add(h[1],h[2],h[3],h[4],h[5],h[6])}return m}function n(c,t){if(null==t){var m=!0;t="linearGradient"==c.type||"radialGradient"==c.type?c.node.getAttribute("gradientTransform"):"pattern"==c.type?c.node.getAttribute("patternTransform"):
c.node.getAttribute("transform");if(!t)return new a.Matrix;t=d(t)}else t=a._.rgTransform.test(t)?J(t).replace(/\.{3}|\u2026/g,c._.transform||aa):d(t),y(t,"array")&&(t=a.path?a.path.toString.call(t):J(t)),c._.transform=t;var b=f(t,c.getBBox(1));if(m)return b;c.matrix=b}function u(c){c=c.node.ownerSVGElement&&x(c.node.ownerSVGElement)||c.node.parentNode&&x(c.node.parentNode)||a.select("svg")||a(0,0);var t=c.select("defs"),t=null==t?!1:t.node;t||(t=r("defs",c.node).node);return t}function p(c){return c.node.ownerSVGElement&&
x(c.node.ownerSVGElement)||a.select("svg")}function b(c,a,m){function b(c){if(null==c)return aa;if(c==+c)return c;v(B,{width:c});try{return B.getBBox().width}catch(a){return 0}}function h(c){if(null==c)return aa;if(c==+c)return c;v(B,{height:c});try{return B.getBBox().height}catch(a){return 0}}function e(b,B){null==a?d[b]=B(c.attr(b)||0):b==a&&(d=B(null==m?c.attr(b)||0:m))}var f=p(c).node,d={},B=f.querySelector(".svg---mgr");B||(B=v("rect"),v(B,{x:-9E9,y:-9E9,width:10,height:10,"class":"svg---mgr",
fill:"none"}),f.appendChild(B));switch(c.type){case "rect":e("rx",b),e("ry",h);case "image":e("width",b),e("height",h);case "text":e("x",b);e("y",h);break;case "circle":e("cx",b);e("cy",h);e("r",b);break;case "ellipse":e("cx",b);e("cy",h);e("rx",b);e("ry",h);break;case "line":e("x1",b);e("x2",b);e("y1",h);e("y2",h);break;case "marker":e("refX",b);e("markerWidth",b);e("refY",h);e("markerHeight",h);break;case "radialGradient":e("fx",b);e("fy",h);break;case "tspan":e("dx",b);e("dy",h);break;default:e(a,
b)}f.removeChild(B);return d}function q(c){y(c,"array")||(c=Array.prototype.slice.call(arguments,0));for(var a=0,b=0,m=this.node;this[a];)delete this[a++];for(a=0;a<c.length;a++)"set"==c[a].type?c[a].forEach(function(c){m.appendChild(c.node)}):m.appendChild(c[a].node);for(var h=m.childNodes,a=0;a<h.length;a++)this[b++]=x(h[a]);return this}function e(c){if(c.snap in E)return E[c.snap];var a=this.id=V(),b;try{b=c.ownerSVGElement}catch(m){}this.node=c;b&&(this.paper=new s(b));this.type=c.tagName;this.anims=
{};this._={transform:[]};c.snap=a;E[a]=this;"g"==this.type&&(this.add=q);if(this.type in{g:1,mask:1,pattern:1})for(var e in s.prototype)s.prototype[h](e)&&(this[e]=s.prototype[e])}function l(c){this.node=c}function r(c,a){var b=v(c);a.appendChild(b);return x(b)}function s(c,a){var b,m,f,d=s.prototype;if(c&&"svg"==c.tagName){if(c.snap in E)return E[c.snap];var l=c.ownerDocument;b=new e(c);m=c.getElementsByTagName("desc")[0];f=c.getElementsByTagName("defs")[0];m||(m=v("desc"),m.appendChild(l.createTextNode("Created with Snap")),
b.node.appendChild(m));f||(f=v("defs"),b.node.appendChild(f));b.defs=f;for(var ca in d)d[h](ca)&&(b[ca]=d[ca]);b.paper=b.root=b}else b=r("svg",G.doc.body),v(b.node,{height:a,version:1.1,width:c,xmlns:la});return b}function x(c){return!c||c instanceof e||c instanceof l?c:c.tagName&&"svg"==c.tagName.toLowerCase()?new s(c):c.tagName&&"object"==c.tagName.toLowerCase()&&"image/svg+xml"==c.type?new s(c.contentDocument.getElementsByTagName("svg")[0]):new e(c)}a.version="0.3.0";a.toString=function(){return"Snap v"+
this.version};a._={};var G={win:N,doc:N.document};a._.glob=G;var h="hasOwnProperty",J=String,K=parseFloat,U=parseInt,I=Math,P=I.max,Q=I.min,Y=I.abs,C=I.PI,aa="",$=Object.prototype.toString,F=/^\s*((#[a-f\d]{6})|(#[a-f\d]{3})|rgba?\(\s*([\d\.]+%?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+%?(?:\s*,\s*[\d\.]+%?)?)\s*\)|hsba?\(\s*([\d\.]+(?:deg|\xb0|%)?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+(?:%?\s*,\s*[\d\.]+)?%?)\s*\)|hsla?\(\s*([\d\.]+(?:deg|\xb0|%)?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+(?:%?\s*,\s*[\d\.]+)?%?)\s*\))\s*$/i;a._.separator=
RegExp("[,\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]+");var S=RegExp("[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*"),X={hs:1,rg:1},W=RegExp("([a-z])[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029,]*((-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*)+)",
"ig"),ma=RegExp("([rstm])[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029,]*((-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*)+)","ig"),Z=RegExp("(-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?)[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*",
"ig"),na=0,ba="S"+(+new Date).toString(36),V=function(){return ba+(na++).toString(36)},m="http://www.w3.org/1999/xlink",la="http://www.w3.org/2000/svg",E={},ca=a.url=function(c){return"url('#"+c+"')"};a._.$=v;a._.id=V;a.format=function(){var c=/\{([^\}]+)\}/g,a=/(?:(?:^|\.)(.+?)(?=\[|\.|$|\()|\[('|")(.+?)\2\])(\(\))?/g,b=function(c,b,m){var h=m;b.replace(a,function(c,a,b,m,t){a=a||m;h&&(a in h&&(h=h[a]),"function"==typeof h&&t&&(h=h()))});return h=(null==h||h==m?c:h)+""};return function(a,m){return J(a).replace(c,
function(c,a){return b(c,a,m)})}}();a._.clone=M;a._.cacher=A;a.rad=z;a.deg=function(c){return 180*c/C%360};a.angle=w;a.is=y;a.snapTo=function(c,a,b){b=y(b,"finite")?b:10;if(y(c,"array"))for(var m=c.length;m--;){if(Y(c[m]-a)<=b)return c[m]}else{c=+c;m=a%c;if(m<b)return a-m;if(m>c-b)return a-m+c}return a};a.getRGB=A(function(c){if(!c||(c=J(c)).indexOf("-")+1)return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka};if("none"==c)return{r:-1,g:-1,b:-1,hex:"none",toString:ka};!X[h](c.toLowerCase().substring(0,
2))&&"#"!=c.charAt()&&(c=T(c));if(!c)return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka};var b,m,e,f,d;if(c=c.match(F)){c[2]&&(e=U(c[2].substring(5),16),m=U(c[2].substring(3,5),16),b=U(c[2].substring(1,3),16));c[3]&&(e=U((d=c[3].charAt(3))+d,16),m=U((d=c[3].charAt(2))+d,16),b=U((d=c[3].charAt(1))+d,16));c[4]&&(d=c[4].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b*=2.55),m=K(d[1]),"%"==d[1].slice(-1)&&(m*=2.55),e=K(d[2]),"%"==d[2].slice(-1)&&(e*=2.55),"rgba"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),
d[3]&&"%"==d[3].slice(-1)&&(f/=100));if(c[5])return d=c[5].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b/=100),m=K(d[1]),"%"==d[1].slice(-1)&&(m/=100),e=K(d[2]),"%"==d[2].slice(-1)&&(e/=100),"deg"!=d[0].slice(-3)&&"\u00b0"!=d[0].slice(-1)||(b/=360),"hsba"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),d[3]&&"%"==d[3].slice(-1)&&(f/=100),a.hsb2rgb(b,m,e,f);if(c[6])return d=c[6].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b/=100),m=K(d[1]),"%"==d[1].slice(-1)&&(m/=100),e=K(d[2]),"%"==d[2].slice(-1)&&(e/=100),
"deg"!=d[0].slice(-3)&&"\u00b0"!=d[0].slice(-1)||(b/=360),"hsla"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),d[3]&&"%"==d[3].slice(-1)&&(f/=100),a.hsl2rgb(b,m,e,f);b=Q(I.round(b),255);m=Q(I.round(m),255);e=Q(I.round(e),255);f=Q(P(f,0),1);c={r:b,g:m,b:e,toString:ka};c.hex="#"+(16777216|e|m<<8|b<<16).toString(16).slice(1);c.opacity=y(f,"finite")?f:1;return c}return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka}},a);a.hsb=A(function(c,b,m){return a.hsb2rgb(c,b,m).hex});a.hsl=A(function(c,b,m){return a.hsl2rgb(c,
b,m).hex});a.rgb=A(function(c,a,b,m){if(y(m,"finite")){var e=I.round;return"rgba("+[e(c),e(a),e(b),+m.toFixed(2)]+")"}return"#"+(16777216|b|a<<8|c<<16).toString(16).slice(1)});var T=function(c){var a=G.doc.getElementsByTagName("head")[0]||G.doc.getElementsByTagName("svg")[0];T=A(function(c){if("red"==c.toLowerCase())return"rgb(255, 0, 0)";a.style.color="rgb(255, 0, 0)";a.style.color=c;c=G.doc.defaultView.getComputedStyle(a,aa).getPropertyValue("color");return"rgb(255, 0, 0)"==c?null:c});return T(c)},
qa=function(){return"hsb("+[this.h,this.s,this.b]+")"},ra=function(){return"hsl("+[this.h,this.s,this.l]+")"},ka=function(){return 1==this.opacity||null==this.opacity?this.hex:"rgba("+[this.r,this.g,this.b,this.opacity]+")"},D=function(c,b,m){null==b&&y(c,"object")&&"r"in c&&"g"in c&&"b"in c&&(m=c.b,b=c.g,c=c.r);null==b&&y(c,string)&&(m=a.getRGB(c),c=m.r,b=m.g,m=m.b);if(1<c||1<b||1<m)c/=255,b/=255,m/=255;return[c,b,m]},oa=function(c,b,m,e){c=I.round(255*c);b=I.round(255*b);m=I.round(255*m);c={r:c,
g:b,b:m,opacity:y(e,"finite")?e:1,hex:a.rgb(c,b,m),toString:ka};y(e,"finite")&&(c.opacity=e);return c};a.color=function(c){var b;y(c,"object")&&"h"in c&&"s"in c&&"b"in c?(b=a.hsb2rgb(c),c.r=b.r,c.g=b.g,c.b=b.b,c.opacity=1,c.hex=b.hex):y(c,"object")&&"h"in c&&"s"in c&&"l"in c?(b=a.hsl2rgb(c),c.r=b.r,c.g=b.g,c.b=b.b,c.opacity=1,c.hex=b.hex):(y(c,"string")&&(c=a.getRGB(c)),y(c,"object")&&"r"in c&&"g"in c&&"b"in c&&!("error"in c)?(b=a.rgb2hsl(c),c.h=b.h,c.s=b.s,c.l=b.l,b=a.rgb2hsb(c),c.v=b.b):(c={hex:"none"},
c.r=c.g=c.b=c.h=c.s=c.v=c.l=-1,c.error=1));c.toString=ka;return c};a.hsb2rgb=function(c,a,b,m){y(c,"object")&&"h"in c&&"s"in c&&"b"in c&&(b=c.b,a=c.s,c=c.h,m=c.o);var e,h,d;c=360*c%360/60;d=b*a;a=d*(1-Y(c%2-1));b=e=h=b-d;c=~~c;b+=[d,a,0,0,a,d][c];e+=[a,d,d,a,0,0][c];h+=[0,0,a,d,d,a][c];return oa(b,e,h,m)};a.hsl2rgb=function(c,a,b,m){y(c,"object")&&"h"in c&&"s"in c&&"l"in c&&(b=c.l,a=c.s,c=c.h);if(1<c||1<a||1<b)c/=360,a/=100,b/=100;var e,h,d;c=360*c%360/60;d=2*a*(0.5>b?b:1-b);a=d*(1-Y(c%2-1));b=e=
h=b-d/2;c=~~c;b+=[d,a,0,0,a,d][c];e+=[a,d,d,a,0,0][c];h+=[0,0,a,d,d,a][c];return oa(b,e,h,m)};a.rgb2hsb=function(c,a,b){b=D(c,a,b);c=b[0];a=b[1];b=b[2];var m,e;m=P(c,a,b);e=m-Q(c,a,b);c=((0==e?0:m==c?(a-b)/e:m==a?(b-c)/e+2:(c-a)/e+4)+360)%6*60/360;return{h:c,s:0==e?0:e/m,b:m,toString:qa}};a.rgb2hsl=function(c,a,b){b=D(c,a,b);c=b[0];a=b[1];b=b[2];var m,e,h;m=P(c,a,b);e=Q(c,a,b);h=m-e;c=((0==h?0:m==c?(a-b)/h:m==a?(b-c)/h+2:(c-a)/h+4)+360)%6*60/360;m=(m+e)/2;return{h:c,s:0==h?0:0.5>m?h/(2*m):h/(2-2*
m),l:m,toString:ra}};a.parsePathString=function(c){if(!c)return null;var b=a.path(c);if(b.arr)return a.path.clone(b.arr);var m={a:7,c:6,o:2,h:1,l:2,m:2,r:4,q:4,s:4,t:2,v:1,u:3,z:0},e=[];y(c,"array")&&y(c[0],"array")&&(e=a.path.clone(c));e.length||J(c).replace(W,function(c,a,b){var h=[];c=a.toLowerCase();b.replace(Z,function(c,a){a&&h.push(+a)});"m"==c&&2<h.length&&(e.push([a].concat(h.splice(0,2))),c="l",a="m"==a?"l":"L");"o"==c&&1==h.length&&e.push([a,h[0] ]);if("r"==c)e.push([a].concat(h));else for(;h.length>=
m[c]&&(e.push([a].concat(h.splice(0,m[c]))),m[c]););});e.toString=a.path.toString;b.arr=a.path.clone(e);return e};var O=a.parseTransformString=function(c){if(!c)return null;var b=[];y(c,"array")&&y(c[0],"array")&&(b=a.path.clone(c));b.length||J(c).replace(ma,function(c,a,m){var e=[];a.toLowerCase();m.replace(Z,function(c,a){a&&e.push(+a)});b.push([a].concat(e))});b.toString=a.path.toString;return b};a._.svgTransform2string=d;a._.rgTransform=RegExp("^[a-z][\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*-?\\.?\\d",
"i");a._.transform2matrix=f;a._unit2px=b;a._.getSomeDefs=u;a._.getSomeSVG=p;a.select=function(c){return x(G.doc.querySelector(c))};a.selectAll=function(c){c=G.doc.querySelectorAll(c);for(var b=(a.set||Array)(),m=0;m<c.length;m++)b.push(x(c[m]));return b};setInterval(function(){for(var c in E)if(E[h](c)){var a=E[c],b=a.node;("svg"!=a.type&&!b.ownerSVGElement||"svg"==a.type&&(!b.parentNode||"ownerSVGElement"in b.parentNode&&!b.ownerSVGElement))&&delete E[c]}},1E4);(function(c){function m(c){function a(c,
b){var m=v(c.node,b);(m=(m=m&&m.match(d))&&m[2])&&"#"==m.charAt()&&(m=m.substring(1))&&(f[m]=(f[m]||[]).concat(function(a){var m={};m[b]=ca(a);v(c.node,m)}))}function b(c){var a=v(c.node,"xlink:href");a&&"#"==a.charAt()&&(a=a.substring(1))&&(f[a]=(f[a]||[]).concat(function(a){c.attr("xlink:href","#"+a)}))}var e=c.selectAll("*"),h,d=/^\s*url\(("|'|)(.*)\1\)\s*$/;c=[];for(var f={},l=0,E=e.length;l<E;l++){h=e[l];a(h,"fill");a(h,"stroke");a(h,"filter");a(h,"mask");a(h,"clip-path");b(h);var t=v(h.node,
"id");t&&(v(h.node,{id:h.id}),c.push({old:t,id:h.id}))}l=0;for(E=c.length;l<E;l++)if(e=f[c[l].old])for(h=0,t=e.length;h<t;h++)e[h](c[l].id)}function e(c,a,b){return function(m){m=m.slice(c,a);1==m.length&&(m=m[0]);return b?b(m):m}}function d(c){return function(){var a=c?"<"+this.type:"",b=this.node.attributes,m=this.node.childNodes;if(c)for(var e=0,h=b.length;e<h;e++)a+=" "+b[e].name+'="'+b[e].value.replace(/"/g,'\\"')+'"';if(m.length){c&&(a+=">");e=0;for(h=m.length;e<h;e++)3==m[e].nodeType?a+=m[e].nodeValue:
1==m[e].nodeType&&(a+=x(m[e]).toString());c&&(a+="</"+this.type+">")}else c&&(a+="/>");return a}}c.attr=function(c,a){if(!c)return this;if(y(c,"string"))if(1<arguments.length){var b={};b[c]=a;c=b}else return k("snap.util.getattr."+c,this).firstDefined();for(var m in c)c[h](m)&&k("snap.util.attr."+m,this,c[m]);return this};c.getBBox=function(c){if(!a.Matrix||!a.path)return this.node.getBBox();var b=this,m=new a.Matrix;if(b.removed)return a._.box();for(;"use"==b.type;)if(c||(m=m.add(b.transform().localMatrix.translate(b.attr("x")||
0,b.attr("y")||0))),b.original)b=b.original;else var e=b.attr("xlink:href"),b=b.original=b.node.ownerDocument.getElementById(e.substring(e.indexOf("#")+1));var e=b._,h=a.path.get[b.type]||a.path.get.deflt;try{if(c)return e.bboxwt=h?a.path.getBBox(b.realPath=h(b)):a._.box(b.node.getBBox()),a._.box(e.bboxwt);b.realPath=h(b);b.matrix=b.transform().localMatrix;e.bbox=a.path.getBBox(a.path.map(b.realPath,m.add(b.matrix)));return a._.box(e.bbox)}catch(d){return a._.box()}};var f=function(){return this.string};
c.transform=function(c){var b=this._;if(null==c){var m=this;c=new a.Matrix(this.node.getCTM());for(var e=n(this),h=[e],d=new a.Matrix,l=e.toTransformString(),b=J(e)==J(this.matrix)?J(b.transform):l;"svg"!=m.type&&(m=m.parent());)h.push(n(m));for(m=h.length;m--;)d.add(h[m]);return{string:b,globalMatrix:c,totalMatrix:d,localMatrix:e,diffMatrix:c.clone().add(e.invert()),global:c.toTransformString(),total:d.toTransformString(),local:l,toString:f}}c instanceof a.Matrix?this.matrix=c:n(this,c);this.node&&
("linearGradient"==this.type||"radialGradient"==this.type?v(this.node,{gradientTransform:this.matrix}):"pattern"==this.type?v(this.node,{patternTransform:this.matrix}):v(this.node,{transform:this.matrix}));return this};c.parent=function(){return x(this.node.parentNode)};c.append=c.add=function(c){if(c){if("set"==c.type){var a=this;c.forEach(function(c){a.add(c)});return this}c=x(c);this.node.appendChild(c.node);c.paper=this.paper}return this};c.appendTo=function(c){c&&(c=x(c),c.append(this));return this};
c.prepend=function(c){if(c){if("set"==c.type){var a=this,b;c.forEach(function(c){b?b.after(c):a.prepend(c);b=c});return this}c=x(c);var m=c.parent();this.node.insertBefore(c.node,this.node.firstChild);this.add&&this.add();c.paper=this.paper;this.parent()&&this.parent().add();m&&m.add()}return this};c.prependTo=function(c){c=x(c);c.prepend(this);return this};c.before=function(c){if("set"==c.type){var a=this;c.forEach(function(c){var b=c.parent();a.node.parentNode.insertBefore(c.node,a.node);b&&b.add()});
this.parent().add();return this}c=x(c);var b=c.parent();this.node.parentNode.insertBefore(c.node,this.node);this.parent()&&this.parent().add();b&&b.add();c.paper=this.paper;return this};c.after=function(c){c=x(c);var a=c.parent();this.node.nextSibling?this.node.parentNode.insertBefore(c.node,this.node.nextSibling):this.node.parentNode.appendChild(c.node);this.parent()&&this.parent().add();a&&a.add();c.paper=this.paper;return this};c.insertBefore=function(c){c=x(c);var a=this.parent();c.node.parentNode.insertBefore(this.node,
c.node);this.paper=c.paper;a&&a.add();c.parent()&&c.parent().add();return this};c.insertAfter=function(c){c=x(c);var a=this.parent();c.node.parentNode.insertBefore(this.node,c.node.nextSibling);this.paper=c.paper;a&&a.add();c.parent()&&c.parent().add();return this};c.remove=function(){var c=this.parent();this.node.parentNode&&this.node.parentNode.removeChild(this.node);delete this.paper;this.removed=!0;c&&c.add();return this};c.select=function(c){return x(this.node.querySelector(c))};c.selectAll=
function(c){c=this.node.querySelectorAll(c);for(var b=(a.set||Array)(),m=0;m<c.length;m++)b.push(x(c[m]));return b};c.asPX=function(c,a){null==a&&(a=this.attr(c));return+b(this,c,a)};c.use=function(){var c,a=this.node.id;a||(a=this.id,v(this.node,{id:a}));c="linearGradient"==this.type||"radialGradient"==this.type||"pattern"==this.type?r(this.type,this.node.parentNode):r("use",this.node.parentNode);v(c.node,{"xlink:href":"#"+a});c.original=this;return c};var l=/\S+/g;c.addClass=function(c){var a=(c||
"").match(l)||[];c=this.node;var b=c.className.baseVal,m=b.match(l)||[],e,h,d;if(a.length){for(e=0;d=a[e++];)h=m.indexOf(d),~h||m.push(d);a=m.join(" ");b!=a&&(c.className.baseVal=a)}return this};c.removeClass=function(c){var a=(c||"").match(l)||[];c=this.node;var b=c.className.baseVal,m=b.match(l)||[],e,h;if(m.length){for(e=0;h=a[e++];)h=m.indexOf(h),~h&&m.splice(h,1);a=m.join(" ");b!=a&&(c.className.baseVal=a)}return this};c.hasClass=function(c){return!!~(this.node.className.baseVal.match(l)||[]).indexOf(c)};
c.toggleClass=function(c,a){if(null!=a)return a?this.addClass(c):this.removeClass(c);var b=(c||"").match(l)||[],m=this.node,e=m.className.baseVal,h=e.match(l)||[],d,f,E;for(d=0;E=b[d++];)f=h.indexOf(E),~f?h.splice(f,1):h.push(E);b=h.join(" ");e!=b&&(m.className.baseVal=b);return this};c.clone=function(){var c=x(this.node.cloneNode(!0));v(c.node,"id")&&v(c.node,{id:c.id});m(c);c.insertAfter(this);return c};c.toDefs=function(){u(this).appendChild(this.node);return this};c.pattern=c.toPattern=function(c,
a,b,m){var e=r("pattern",u(this));null==c&&(c=this.getBBox());y(c,"object")&&"x"in c&&(a=c.y,b=c.width,m=c.height,c=c.x);v(e.node,{x:c,y:a,width:b,height:m,patternUnits:"userSpaceOnUse",id:e.id,viewBox:[c,a,b,m].join(" ")});e.node.appendChild(this.node);return e};c.marker=function(c,a,b,m,e,h){var d=r("marker",u(this));null==c&&(c=this.getBBox());y(c,"object")&&"x"in c&&(a=c.y,b=c.width,m=c.height,e=c.refX||c.cx,h=c.refY||c.cy,c=c.x);v(d.node,{viewBox:[c,a,b,m].join(" "),markerWidth:b,markerHeight:m,
orient:"auto",refX:e||0,refY:h||0,id:d.id});d.node.appendChild(this.node);return d};var E=function(c,a,b,m){"function"!=typeof b||b.length||(m=b,b=L.linear);this.attr=c;this.dur=a;b&&(this.easing=b);m&&(this.callback=m)};a._.Animation=E;a.animation=function(c,a,b,m){return new E(c,a,b,m)};c.inAnim=function(){var c=[],a;for(a in this.anims)this.anims[h](a)&&function(a){c.push({anim:new E(a._attrs,a.dur,a.easing,a._callback),mina:a,curStatus:a.status(),status:function(c){return a.status(c)},stop:function(){a.stop()}})}(this.anims[a]);
return c};a.animate=function(c,a,b,m,e,h){"function"!=typeof e||e.length||(h=e,e=L.linear);var d=L.time();c=L(c,a,d,d+m,L.time,b,e);h&&k.once("mina.finish."+c.id,h);return c};c.stop=function(){for(var c=this.inAnim(),a=0,b=c.length;a<b;a++)c[a].stop();return this};c.animate=function(c,a,b,m){"function"!=typeof b||b.length||(m=b,b=L.linear);c instanceof E&&(m=c.callback,b=c.easing,a=b.dur,c=c.attr);var d=[],f=[],l={},t,ca,n,T=this,q;for(q in c)if(c[h](q)){T.equal?(n=T.equal(q,J(c[q])),t=n.from,ca=
n.to,n=n.f):(t=+T.attr(q),ca=+c[q]);var la=y(t,"array")?t.length:1;l[q]=e(d.length,d.length+la,n);d=d.concat(t);f=f.concat(ca)}t=L.time();var p=L(d,f,t,t+a,L.time,function(c){var a={},b;for(b in l)l[h](b)&&(a[b]=l[b](c));T.attr(a)},b);T.anims[p.id]=p;p._attrs=c;p._callback=m;k("snap.animcreated."+T.id,p);k.once("mina.finish."+p.id,function(){delete T.anims[p.id];m&&m.call(T)});k.once("mina.stop."+p.id,function(){delete T.anims[p.id]});return T};var T={};c.data=function(c,b){var m=T[this.id]=T[this.id]||
{};if(0==arguments.length)return k("snap.data.get."+this.id,this,m,null),m;if(1==arguments.length){if(a.is(c,"object")){for(var e in c)c[h](e)&&this.data(e,c[e]);return this}k("snap.data.get."+this.id,this,m[c],c);return m[c]}m[c]=b;k("snap.data.set."+this.id,this,b,c);return this};c.removeData=function(c){null==c?T[this.id]={}:T[this.id]&&delete T[this.id][c];return this};c.outerSVG=c.toString=d(1);c.innerSVG=d()})(e.prototype);a.parse=function(c){var a=G.doc.createDocumentFragment(),b=!0,m=G.doc.createElement("div");
c=J(c);c.match(/^\s*<\s*svg(?:\s|>)/)||(c="<svg>"+c+"</svg>",b=!1);m.innerHTML=c;if(c=m.getElementsByTagName("svg")[0])if(b)a=c;else for(;c.firstChild;)a.appendChild(c.firstChild);m.innerHTML=aa;return new l(a)};l.prototype.select=e.prototype.select;l.prototype.selectAll=e.prototype.selectAll;a.fragment=function(){for(var c=Array.prototype.slice.call(arguments,0),b=G.doc.createDocumentFragment(),m=0,e=c.length;m<e;m++){var h=c[m];h.node&&h.node.nodeType&&b.appendChild(h.node);h.nodeType&&b.appendChild(h);
"string"==typeof h&&b.appendChild(a.parse(h).node)}return new l(b)};a._.make=r;a._.wrap=x;s.prototype.el=function(c,a){var b=r(c,this.node);a&&b.attr(a);return b};k.on("snap.util.getattr",function(){var c=k.nt(),c=c.substring(c.lastIndexOf(".")+1),a=c.replace(/[A-Z]/g,function(c){return"-"+c.toLowerCase()});return pa[h](a)?this.node.ownerDocument.defaultView.getComputedStyle(this.node,null).getPropertyValue(a):v(this.node,c)});var pa={"alignment-baseline":0,"baseline-shift":0,clip:0,"clip-path":0,
"clip-rule":0,color:0,"color-interpolation":0,"color-interpolation-filters":0,"color-profile":0,"color-rendering":0,cursor:0,direction:0,display:0,"dominant-baseline":0,"enable-background":0,fill:0,"fill-opacity":0,"fill-rule":0,filter:0,"flood-color":0,"flood-opacity":0,font:0,"font-family":0,"font-size":0,"font-size-adjust":0,"font-stretch":0,"font-style":0,"font-variant":0,"font-weight":0,"glyph-orientation-horizontal":0,"glyph-orientation-vertical":0,"image-rendering":0,kerning:0,"letter-spacing":0,
"lighting-color":0,marker:0,"marker-end":0,"marker-mid":0,"marker-start":0,mask:0,opacity:0,overflow:0,"pointer-events":0,"shape-rendering":0,"stop-color":0,"stop-opacity":0,stroke:0,"stroke-dasharray":0,"stroke-dashoffset":0,"stroke-linecap":0,"stroke-linejoin":0,"stroke-miterlimit":0,"stroke-opacity":0,"stroke-width":0,"text-anchor":0,"text-decoration":0,"text-rendering":0,"unicode-bidi":0,visibility:0,"word-spacing":0,"writing-mode":0};k.on("snap.util.attr",function(c){var a=k.nt(),b={},a=a.substring(a.lastIndexOf(".")+
1);b[a]=c;var m=a.replace(/-(\w)/gi,function(c,a){return a.toUpperCase()}),a=a.replace(/[A-Z]/g,function(c){return"-"+c.toLowerCase()});pa[h](a)?this.node.style[m]=null==c?aa:c:v(this.node,b)});a.ajax=function(c,a,b,m){var e=new XMLHttpRequest,h=V();if(e){if(y(a,"function"))m=b,b=a,a=null;else if(y(a,"object")){var d=[],f;for(f in a)a.hasOwnProperty(f)&&d.push(encodeURIComponent(f)+"="+encodeURIComponent(a[f]));a=d.join("&")}e.open(a?"POST":"GET",c,!0);a&&(e.setRequestHeader("X-Requested-With","XMLHttpRequest"),
e.setRequestHeader("Content-type","application/x-www-form-urlencoded"));b&&(k.once("snap.ajax."+h+".0",b),k.once("snap.ajax."+h+".200",b),k.once("snap.ajax."+h+".304",b));e.onreadystatechange=function(){4==e.readyState&&k("snap.ajax."+h+"."+e.status,m,e)};if(4==e.readyState)return e;e.send(a);return e}};a.load=function(c,b,m){a.ajax(c,function(c){c=a.parse(c.responseText);m?b.call(m,c):b(c)})};a.getElementByPoint=function(c,a){var b,m,e=G.doc.elementFromPoint(c,a);if(G.win.opera&&"svg"==e.tagName){b=
e;m=b.getBoundingClientRect();b=b.ownerDocument;var h=b.body,d=b.documentElement;b=m.top+(g.win.pageYOffset||d.scrollTop||h.scrollTop)-(d.clientTop||h.clientTop||0);m=m.left+(g.win.pageXOffset||d.scrollLeft||h.scrollLeft)-(d.clientLeft||h.clientLeft||0);h=e.createSVGRect();h.x=c-m;h.y=a-b;h.width=h.height=1;b=e.getIntersectionList(h,null);b.length&&(e=b[b.length-1])}return e?x(e):null};a.plugin=function(c){c(a,e,s,G,l)};return G.win.Snap=a}();C.plugin(function(a,k,y,M,A){function w(a,d,f,b,q,e){null==
d&&"[object SVGMatrix]"==z.call(a)?(this.a=a.a,this.b=a.b,this.c=a.c,this.d=a.d,this.e=a.e,this.f=a.f):null!=a?(this.a=+a,this.b=+d,this.c=+f,this.d=+b,this.e=+q,this.f=+e):(this.a=1,this.c=this.b=0,this.d=1,this.f=this.e=0)}var z=Object.prototype.toString,d=String,f=Math;(function(n){function k(a){return a[0]*a[0]+a[1]*a[1]}function p(a){var d=f.sqrt(k(a));a[0]&&(a[0]/=d);a[1]&&(a[1]/=d)}n.add=function(a,d,e,f,n,p){var k=[[],[],[] ],u=[[this.a,this.c,this.e],[this.b,this.d,this.f],[0,0,1] ];d=[[a,
e,n],[d,f,p],[0,0,1] ];a&&a instanceof w&&(d=[[a.a,a.c,a.e],[a.b,a.d,a.f],[0,0,1] ]);for(a=0;3>a;a++)for(e=0;3>e;e++){for(f=n=0;3>f;f++)n+=u[a][f]*d[f][e];k[a][e]=n}this.a=k[0][0];this.b=k[1][0];this.c=k[0][1];this.d=k[1][1];this.e=k[0][2];this.f=k[1][2];return this};n.invert=function(){var a=this.a*this.d-this.b*this.c;return new w(this.d/a,-this.b/a,-this.c/a,this.a/a,(this.c*this.f-this.d*this.e)/a,(this.b*this.e-this.a*this.f)/a)};n.clone=function(){return new w(this.a,this.b,this.c,this.d,this.e,
this.f)};n.translate=function(a,d){return this.add(1,0,0,1,a,d)};n.scale=function(a,d,e,f){null==d&&(d=a);(e||f)&&this.add(1,0,0,1,e,f);this.add(a,0,0,d,0,0);(e||f)&&this.add(1,0,0,1,-e,-f);return this};n.rotate=function(b,d,e){b=a.rad(b);d=d||0;e=e||0;var l=+f.cos(b).toFixed(9);b=+f.sin(b).toFixed(9);this.add(l,b,-b,l,d,e);return this.add(1,0,0,1,-d,-e)};n.x=function(a,d){return a*this.a+d*this.c+this.e};n.y=function(a,d){return a*this.b+d*this.d+this.f};n.get=function(a){return+this[d.fromCharCode(97+
a)].toFixed(4)};n.toString=function(){return"matrix("+[this.get(0),this.get(1),this.get(2),this.get(3),this.get(4),this.get(5)].join()+")"};n.offset=function(){return[this.e.toFixed(4),this.f.toFixed(4)]};n.determinant=function(){return this.a*this.d-this.b*this.c};n.split=function(){var b={};b.dx=this.e;b.dy=this.f;var d=[[this.a,this.c],[this.b,this.d] ];b.scalex=f.sqrt(k(d[0]));p(d[0]);b.shear=d[0][0]*d[1][0]+d[0][1]*d[1][1];d[1]=[d[1][0]-d[0][0]*b.shear,d[1][1]-d[0][1]*b.shear];b.scaley=f.sqrt(k(d[1]));
p(d[1]);b.shear/=b.scaley;0>this.determinant()&&(b.scalex=-b.scalex);var e=-d[0][1],d=d[1][1];0>d?(b.rotate=a.deg(f.acos(d)),0>e&&(b.rotate=360-b.rotate)):b.rotate=a.deg(f.asin(e));b.isSimple=!+b.shear.toFixed(9)&&(b.scalex.toFixed(9)==b.scaley.toFixed(9)||!b.rotate);b.isSuperSimple=!+b.shear.toFixed(9)&&b.scalex.toFixed(9)==b.scaley.toFixed(9)&&!b.rotate;b.noRotation=!+b.shear.toFixed(9)&&!b.rotate;return b};n.toTransformString=function(a){a=a||this.split();if(+a.shear.toFixed(9))return"m"+[this.get(0),
this.get(1),this.get(2),this.get(3),this.get(4),this.get(5)];a.scalex=+a.scalex.toFixed(4);a.scaley=+a.scaley.toFixed(4);a.rotate=+a.rotate.toFixed(4);return(a.dx||a.dy?"t"+[+a.dx.toFixed(4),+a.dy.toFixed(4)]:"")+(1!=a.scalex||1!=a.scaley?"s"+[a.scalex,a.scaley,0,0]:"")+(a.rotate?"r"+[+a.rotate.toFixed(4),0,0]:"")}})(w.prototype);a.Matrix=w;a.matrix=function(a,d,f,b,k,e){return new w(a,d,f,b,k,e)}});C.plugin(function(a,v,y,M,A){function w(h){return function(d){k.stop();d instanceof A&&1==d.node.childNodes.length&&
("radialGradient"==d.node.firstChild.tagName||"linearGradient"==d.node.firstChild.tagName||"pattern"==d.node.firstChild.tagName)&&(d=d.node.firstChild,b(this).appendChild(d),d=u(d));if(d instanceof v)if("radialGradient"==d.type||"linearGradient"==d.type||"pattern"==d.type){d.node.id||e(d.node,{id:d.id});var f=l(d.node.id)}else f=d.attr(h);else f=a.color(d),f.error?(f=a(b(this).ownerSVGElement).gradient(d))?(f.node.id||e(f.node,{id:f.id}),f=l(f.node.id)):f=d:f=r(f);d={};d[h]=f;e(this.node,d);this.node.style[h]=
x}}function z(a){k.stop();a==+a&&(a+="px");this.node.style.fontSize=a}function d(a){var b=[];a=a.childNodes;for(var e=0,f=a.length;e<f;e++){var l=a[e];3==l.nodeType&&b.push(l.nodeValue);"tspan"==l.tagName&&(1==l.childNodes.length&&3==l.firstChild.nodeType?b.push(l.firstChild.nodeValue):b.push(d(l)))}return b}function f(){k.stop();return this.node.style.fontSize}var n=a._.make,u=a._.wrap,p=a.is,b=a._.getSomeDefs,q=/^url\(#?([^)]+)\)$/,e=a._.$,l=a.url,r=String,s=a._.separator,x="";k.on("snap.util.attr.mask",
function(a){if(a instanceof v||a instanceof A){k.stop();a instanceof A&&1==a.node.childNodes.length&&(a=a.node.firstChild,b(this).appendChild(a),a=u(a));if("mask"==a.type)var d=a;else d=n("mask",b(this)),d.node.appendChild(a.node);!d.node.id&&e(d.node,{id:d.id});e(this.node,{mask:l(d.id)})}});(function(a){k.on("snap.util.attr.clip",a);k.on("snap.util.attr.clip-path",a);k.on("snap.util.attr.clipPath",a)})(function(a){if(a instanceof v||a instanceof A){k.stop();if("clipPath"==a.type)var d=a;else d=
n("clipPath",b(this)),d.node.appendChild(a.node),!d.node.id&&e(d.node,{id:d.id});e(this.node,{"clip-path":l(d.id)})}});k.on("snap.util.attr.fill",w("fill"));k.on("snap.util.attr.stroke",w("stroke"));var G=/^([lr])(?:\(([^)]*)\))?(.*)$/i;k.on("snap.util.grad.parse",function(a){a=r(a);var b=a.match(G);if(!b)return null;a=b[1];var e=b[2],b=b[3],e=e.split(/\s*,\s*/).map(function(a){return+a==a?+a:a});1==e.length&&0==e[0]&&(e=[]);b=b.split("-");b=b.map(function(a){a=a.split(":");var b={color:a[0]};a[1]&&
(b.offset=parseFloat(a[1]));return b});return{type:a,params:e,stops:b}});k.on("snap.util.attr.d",function(b){k.stop();p(b,"array")&&p(b[0],"array")&&(b=a.path.toString.call(b));b=r(b);b.match(/[ruo]/i)&&(b=a.path.toAbsolute(b));e(this.node,{d:b})})(-1);k.on("snap.util.attr.#text",function(a){k.stop();a=r(a);for(a=M.doc.createTextNode(a);this.node.firstChild;)this.node.removeChild(this.node.firstChild);this.node.appendChild(a)})(-1);k.on("snap.util.attr.path",function(a){k.stop();this.attr({d:a})})(-1);
k.on("snap.util.attr.class",function(a){k.stop();this.node.className.baseVal=a})(-1);k.on("snap.util.attr.viewBox",function(a){a=p(a,"object")&&"x"in a?[a.x,a.y,a.width,a.height].join(" "):p(a,"array")?a.join(" "):a;e(this.node,{viewBox:a});k.stop()})(-1);k.on("snap.util.attr.transform",function(a){this.transform(a);k.stop()})(-1);k.on("snap.util.attr.r",function(a){"rect"==this.type&&(k.stop(),e(this.node,{rx:a,ry:a}))})(-1);k.on("snap.util.attr.textpath",function(a){k.stop();if("text"==this.type){var d,
f;if(!a&&this.textPath){for(a=this.textPath;a.node.firstChild;)this.node.appendChild(a.node.firstChild);a.remove();delete this.textPath}else if(p(a,"string")?(d=b(this),a=u(d.parentNode).path(a),d.appendChild(a.node),d=a.id,a.attr({id:d})):(a=u(a),a instanceof v&&(d=a.attr("id"),d||(d=a.id,a.attr({id:d})))),d)if(a=this.textPath,f=this.node,a)a.attr({"xlink:href":"#"+d});else{for(a=e("textPath",{"xlink:href":"#"+d});f.firstChild;)a.appendChild(f.firstChild);f.appendChild(a);this.textPath=u(a)}}})(-1);
k.on("snap.util.attr.text",function(a){if("text"==this.type){for(var b=this.node,d=function(a){var b=e("tspan");if(p(a,"array"))for(var f=0;f<a.length;f++)b.appendChild(d(a[f]));else b.appendChild(M.doc.createTextNode(a));b.normalize&&b.normalize();return b};b.firstChild;)b.removeChild(b.firstChild);for(a=d(a);a.firstChild;)b.appendChild(a.firstChild)}k.stop()})(-1);k.on("snap.util.attr.fontSize",z)(-1);k.on("snap.util.attr.font-size",z)(-1);k.on("snap.util.getattr.transform",function(){k.stop();
return this.transform()})(-1);k.on("snap.util.getattr.textpath",function(){k.stop();return this.textPath})(-1);(function(){function b(d){return function(){k.stop();var b=M.doc.defaultView.getComputedStyle(this.node,null).getPropertyValue("marker-"+d);return"none"==b?b:a(M.doc.getElementById(b.match(q)[1]))}}function d(a){return function(b){k.stop();var d="marker"+a.charAt(0).toUpperCase()+a.substring(1);if(""==b||!b)this.node.style[d]="none";else if("marker"==b.type){var f=b.node.id;f||e(b.node,{id:b.id});
this.node.style[d]=l(f)}}}k.on("snap.util.getattr.marker-end",b("end"))(-1);k.on("snap.util.getattr.markerEnd",b("end"))(-1);k.on("snap.util.getattr.marker-start",b("start"))(-1);k.on("snap.util.getattr.markerStart",b("start"))(-1);k.on("snap.util.getattr.marker-mid",b("mid"))(-1);k.on("snap.util.getattr.markerMid",b("mid"))(-1);k.on("snap.util.attr.marker-end",d("end"))(-1);k.on("snap.util.attr.markerEnd",d("end"))(-1);k.on("snap.util.attr.marker-start",d("start"))(-1);k.on("snap.util.attr.markerStart",
d("start"))(-1);k.on("snap.util.attr.marker-mid",d("mid"))(-1);k.on("snap.util.attr.markerMid",d("mid"))(-1)})();k.on("snap.util.getattr.r",function(){if("rect"==this.type&&e(this.node,"rx")==e(this.node,"ry"))return k.stop(),e(this.node,"rx")})(-1);k.on("snap.util.getattr.text",function(){if("text"==this.type||"tspan"==this.type){k.stop();var a=d(this.node);return 1==a.length?a[0]:a}})(-1);k.on("snap.util.getattr.#text",function(){return this.node.textContent})(-1);k.on("snap.util.getattr.viewBox",
function(){k.stop();var b=e(this.node,"viewBox");if(b)return b=b.split(s),a._.box(+b[0],+b[1],+b[2],+b[3])})(-1);k.on("snap.util.getattr.points",function(){var a=e(this.node,"points");k.stop();if(a)return a.split(s)})(-1);k.on("snap.util.getattr.path",function(){var a=e(this.node,"d");k.stop();return a})(-1);k.on("snap.util.getattr.class",function(){return this.node.className.baseVal})(-1);k.on("snap.util.getattr.fontSize",f)(-1);k.on("snap.util.getattr.font-size",f)(-1)});C.plugin(function(a,v,y,
M,A){function w(a){return a}function z(a){return function(b){return+b.toFixed(3)+a}}var d={"+":function(a,b){return a+b},"-":function(a,b){return a-b},"/":function(a,b){return a/b},"*":function(a,b){return a*b}},f=String,n=/[a-z]+$/i,u=/^\s*([+\-\/*])\s*=\s*([\d.eE+\-]+)\s*([^\d\s]+)?\s*$/;k.on("snap.util.attr",function(a){if(a=f(a).match(u)){var b=k.nt(),b=b.substring(b.lastIndexOf(".")+1),q=this.attr(b),e={};k.stop();var l=a[3]||"",r=q.match(n),s=d[a[1] ];r&&r==l?a=s(parseFloat(q),+a[2]):(q=this.asPX(b),
a=s(this.asPX(b),this.asPX(b,a[2]+l)));isNaN(q)||isNaN(a)||(e[b]=a,this.attr(e))}})(-10);k.on("snap.util.equal",function(a,b){var q=f(this.attr(a)||""),e=f(b).match(u);if(e){k.stop();var l=e[3]||"",r=q.match(n),s=d[e[1] ];if(r&&r==l)return{from:parseFloat(q),to:s(parseFloat(q),+e[2]),f:z(r)};q=this.asPX(a);return{from:q,to:s(q,this.asPX(a,e[2]+l)),f:w}}})(-10)});C.plugin(function(a,v,y,M,A){var w=y.prototype,z=a.is;w.rect=function(a,d,k,p,b,q){var e;null==q&&(q=b);z(a,"object")&&"[object Object]"==
a?e=a:null!=a&&(e={x:a,y:d,width:k,height:p},null!=b&&(e.rx=b,e.ry=q));return this.el("rect",e)};w.circle=function(a,d,k){var p;z(a,"object")&&"[object Object]"==a?p=a:null!=a&&(p={cx:a,cy:d,r:k});return this.el("circle",p)};var d=function(){function a(){this.parentNode.removeChild(this)}return function(d,k){var p=M.doc.createElement("img"),b=M.doc.body;p.style.cssText="position:absolute;left:-9999em;top:-9999em";p.onload=function(){k.call(p);p.onload=p.onerror=null;b.removeChild(p)};p.onerror=a;
b.appendChild(p);p.src=d}}();w.image=function(f,n,k,p,b){var q=this.el("image");if(z(f,"object")&&"src"in f)q.attr(f);else if(null!=f){var e={"xlink:href":f,preserveAspectRatio:"none"};null!=n&&null!=k&&(e.x=n,e.y=k);null!=p&&null!=b?(e.width=p,e.height=b):d(f,function(){a._.$(q.node,{width:this.offsetWidth,height:this.offsetHeight})});a._.$(q.node,e)}return q};w.ellipse=function(a,d,k,p){var b;z(a,"object")&&"[object Object]"==a?b=a:null!=a&&(b={cx:a,cy:d,rx:k,ry:p});return this.el("ellipse",b)};
w.path=function(a){var d;z(a,"object")&&!z(a,"array")?d=a:a&&(d={d:a});return this.el("path",d)};w.group=w.g=function(a){var d=this.el("g");1==arguments.length&&a&&!a.type?d.attr(a):arguments.length&&d.add(Array.prototype.slice.call(arguments,0));return d};w.svg=function(a,d,k,p,b,q,e,l){var r={};z(a,"object")&&null==d?r=a:(null!=a&&(r.x=a),null!=d&&(r.y=d),null!=k&&(r.width=k),null!=p&&(r.height=p),null!=b&&null!=q&&null!=e&&null!=l&&(r.viewBox=[b,q,e,l]));return this.el("svg",r)};w.mask=function(a){var d=
this.el("mask");1==arguments.length&&a&&!a.type?d.attr(a):arguments.length&&d.add(Array.prototype.slice.call(arguments,0));return d};w.ptrn=function(a,d,k,p,b,q,e,l){if(z(a,"object"))var r=a;else arguments.length?(r={},null!=a&&(r.x=a),null!=d&&(r.y=d),null!=k&&(r.width=k),null!=p&&(r.height=p),null!=b&&null!=q&&null!=e&&null!=l&&(r.viewBox=[b,q,e,l])):r={patternUnits:"userSpaceOnUse"};return this.el("pattern",r)};w.use=function(a){return null!=a?(make("use",this.node),a instanceof v&&(a.attr("id")||
a.attr({id:ID()}),a=a.attr("id")),this.el("use",{"xlink:href":a})):v.prototype.use.call(this)};w.text=function(a,d,k){var p={};z(a,"object")?p=a:null!=a&&(p={x:a,y:d,text:k||""});return this.el("text",p)};w.line=function(a,d,k,p){var b={};z(a,"object")?b=a:null!=a&&(b={x1:a,x2:k,y1:d,y2:p});return this.el("line",b)};w.polyline=function(a){1<arguments.length&&(a=Array.prototype.slice.call(arguments,0));var d={};z(a,"object")&&!z(a,"array")?d=a:null!=a&&(d={points:a});return this.el("polyline",d)};
w.polygon=function(a){1<arguments.length&&(a=Array.prototype.slice.call(arguments,0));var d={};z(a,"object")&&!z(a,"array")?d=a:null!=a&&(d={points:a});return this.el("polygon",d)};(function(){function d(){return this.selectAll("stop")}function n(b,d){var f=e("stop"),k={offset:+d+"%"};b=a.color(b);k["stop-color"]=b.hex;1>b.opacity&&(k["stop-opacity"]=b.opacity);e(f,k);this.node.appendChild(f);return this}function u(){if("linearGradient"==this.type){var b=e(this.node,"x1")||0,d=e(this.node,"x2")||
1,f=e(this.node,"y1")||0,k=e(this.node,"y2")||0;return a._.box(b,f,math.abs(d-b),math.abs(k-f))}b=this.node.r||0;return a._.box((this.node.cx||0.5)-b,(this.node.cy||0.5)-b,2*b,2*b)}function p(a,d){function f(a,b){for(var d=(b-u)/(a-w),e=w;e<a;e++)h[e].offset=+(+u+d*(e-w)).toFixed(2);w=a;u=b}var n=k("snap.util.grad.parse",null,d).firstDefined(),p;if(!n)return null;n.params.unshift(a);p="l"==n.type.toLowerCase()?b.apply(0,n.params):q.apply(0,n.params);n.type!=n.type.toLowerCase()&&e(p.node,{gradientUnits:"userSpaceOnUse"});
var h=n.stops,n=h.length,u=0,w=0;n--;for(var v=0;v<n;v++)"offset"in h[v]&&f(v,h[v].offset);h[n].offset=h[n].offset||100;f(n,h[n].offset);for(v=0;v<=n;v++){var y=h[v];p.addStop(y.color,y.offset)}return p}function b(b,k,p,q,w){b=a._.make("linearGradient",b);b.stops=d;b.addStop=n;b.getBBox=u;null!=k&&e(b.node,{x1:k,y1:p,x2:q,y2:w});return b}function q(b,k,p,q,w,h){b=a._.make("radialGradient",b);b.stops=d;b.addStop=n;b.getBBox=u;null!=k&&e(b.node,{cx:k,cy:p,r:q});null!=w&&null!=h&&e(b.node,{fx:w,fy:h});
return b}var e=a._.$;w.gradient=function(a){return p(this.defs,a)};w.gradientLinear=function(a,d,e,f){return b(this.defs,a,d,e,f)};w.gradientRadial=function(a,b,d,e,f){return q(this.defs,a,b,d,e,f)};w.toString=function(){var b=this.node.ownerDocument,d=b.createDocumentFragment(),b=b.createElement("div"),e=this.node.cloneNode(!0);d.appendChild(b);b.appendChild(e);a._.$(e,{xmlns:"http://www.w3.org/2000/svg"});b=b.innerHTML;d.removeChild(d.firstChild);return b};w.clear=function(){for(var a=this.node.firstChild,
b;a;)b=a.nextSibling,"defs"!=a.tagName?a.parentNode.removeChild(a):w.clear.call({node:a}),a=b}})()});C.plugin(function(a,k,y,M){function A(a){var b=A.ps=A.ps||{};b[a]?b[a].sleep=100:b[a]={sleep:100};setTimeout(function(){for(var d in b)b[L](d)&&d!=a&&(b[d].sleep--,!b[d].sleep&&delete b[d])});return b[a]}function w(a,b,d,e){null==a&&(a=b=d=e=0);null==b&&(b=a.y,d=a.width,e=a.height,a=a.x);return{x:a,y:b,width:d,w:d,height:e,h:e,x2:a+d,y2:b+e,cx:a+d/2,cy:b+e/2,r1:F.min(d,e)/2,r2:F.max(d,e)/2,r0:F.sqrt(d*
d+e*e)/2,path:s(a,b,d,e),vb:[a,b,d,e].join(" ")}}function z(){return this.join(",").replace(N,"$1")}function d(a){a=C(a);a.toString=z;return a}function f(a,b,d,h,f,k,l,n,p){if(null==p)return e(a,b,d,h,f,k,l,n);if(0>p||e(a,b,d,h,f,k,l,n)<p)p=void 0;else{var q=0.5,O=1-q,s;for(s=e(a,b,d,h,f,k,l,n,O);0.01<Z(s-p);)q/=2,O+=(s<p?1:-1)*q,s=e(a,b,d,h,f,k,l,n,O);p=O}return u(a,b,d,h,f,k,l,n,p)}function n(b,d){function e(a){return+(+a).toFixed(3)}return a._.cacher(function(a,h,l){a instanceof k&&(a=a.attr("d"));
a=I(a);for(var n,p,D,q,O="",s={},c=0,t=0,r=a.length;t<r;t++){D=a[t];if("M"==D[0])n=+D[1],p=+D[2];else{q=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6]);if(c+q>h){if(d&&!s.start){n=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6],h-c);O+=["C"+e(n.start.x),e(n.start.y),e(n.m.x),e(n.m.y),e(n.x),e(n.y)];if(l)return O;s.start=O;O=["M"+e(n.x),e(n.y)+"C"+e(n.n.x),e(n.n.y),e(n.end.x),e(n.end.y),e(D[5]),e(D[6])].join();c+=q;n=+D[5];p=+D[6];continue}if(!b&&!d)return n=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6],h-c)}c+=q;n=+D[5];p=+D[6]}O+=
D.shift()+D}s.end=O;return n=b?c:d?s:u(n,p,D[0],D[1],D[2],D[3],D[4],D[5],1)},null,a._.clone)}function u(a,b,d,e,h,f,k,l,n){var p=1-n,q=ma(p,3),s=ma(p,2),c=n*n,t=c*n,r=q*a+3*s*n*d+3*p*n*n*h+t*k,q=q*b+3*s*n*e+3*p*n*n*f+t*l,s=a+2*n*(d-a)+c*(h-2*d+a),t=b+2*n*(e-b)+c*(f-2*e+b),x=d+2*n*(h-d)+c*(k-2*h+d),c=e+2*n*(f-e)+c*(l-2*f+e);a=p*a+n*d;b=p*b+n*e;h=p*h+n*k;f=p*f+n*l;l=90-180*F.atan2(s-x,t-c)/S;return{x:r,y:q,m:{x:s,y:t},n:{x:x,y:c},start:{x:a,y:b},end:{x:h,y:f},alpha:l}}function p(b,d,e,h,f,n,k,l){a.is(b,
"array")||(b=[b,d,e,h,f,n,k,l]);b=U.apply(null,b);return w(b.min.x,b.min.y,b.max.x-b.min.x,b.max.y-b.min.y)}function b(a,b,d){return b>=a.x&&b<=a.x+a.width&&d>=a.y&&d<=a.y+a.height}function q(a,d){a=w(a);d=w(d);return b(d,a.x,a.y)||b(d,a.x2,a.y)||b(d,a.x,a.y2)||b(d,a.x2,a.y2)||b(a,d.x,d.y)||b(a,d.x2,d.y)||b(a,d.x,d.y2)||b(a,d.x2,d.y2)||(a.x<d.x2&&a.x>d.x||d.x<a.x2&&d.x>a.x)&&(a.y<d.y2&&a.y>d.y||d.y<a.y2&&d.y>a.y)}function e(a,b,d,e,h,f,n,k,l){null==l&&(l=1);l=(1<l?1:0>l?0:l)/2;for(var p=[-0.1252,
0.1252,-0.3678,0.3678,-0.5873,0.5873,-0.7699,0.7699,-0.9041,0.9041,-0.9816,0.9816],q=[0.2491,0.2491,0.2335,0.2335,0.2032,0.2032,0.1601,0.1601,0.1069,0.1069,0.0472,0.0472],s=0,c=0;12>c;c++)var t=l*p[c]+l,r=t*(t*(-3*a+9*d-9*h+3*n)+6*a-12*d+6*h)-3*a+3*d,t=t*(t*(-3*b+9*e-9*f+3*k)+6*b-12*e+6*f)-3*b+3*e,s=s+q[c]*F.sqrt(r*r+t*t);return l*s}function l(a,b,d){a=I(a);b=I(b);for(var h,f,l,n,k,s,r,O,x,c,t=d?0:[],w=0,v=a.length;w<v;w++)if(x=a[w],"M"==x[0])h=k=x[1],f=s=x[2];else{"C"==x[0]?(x=[h,f].concat(x.slice(1)),
h=x[6],f=x[7]):(x=[h,f,h,f,k,s,k,s],h=k,f=s);for(var G=0,y=b.length;G<y;G++)if(c=b[G],"M"==c[0])l=r=c[1],n=O=c[2];else{"C"==c[0]?(c=[l,n].concat(c.slice(1)),l=c[6],n=c[7]):(c=[l,n,l,n,r,O,r,O],l=r,n=O);var z;var K=x,B=c;z=d;var H=p(K),J=p(B);if(q(H,J)){for(var H=e.apply(0,K),J=e.apply(0,B),H=~~(H/8),J=~~(J/8),U=[],A=[],F={},M=z?0:[],P=0;P<H+1;P++){var C=u.apply(0,K.concat(P/H));U.push({x:C.x,y:C.y,t:P/H})}for(P=0;P<J+1;P++)C=u.apply(0,B.concat(P/J)),A.push({x:C.x,y:C.y,t:P/J});for(P=0;P<H;P++)for(K=
0;K<J;K++){var Q=U[P],L=U[P+1],B=A[K],C=A[K+1],N=0.001>Z(L.x-Q.x)?"y":"x",S=0.001>Z(C.x-B.x)?"y":"x",R;R=Q.x;var Y=Q.y,V=L.x,ea=L.y,fa=B.x,ga=B.y,ha=C.x,ia=C.y;if(W(R,V)<X(fa,ha)||X(R,V)>W(fa,ha)||W(Y,ea)<X(ga,ia)||X(Y,ea)>W(ga,ia))R=void 0;else{var $=(R*ea-Y*V)*(fa-ha)-(R-V)*(fa*ia-ga*ha),aa=(R*ea-Y*V)*(ga-ia)-(Y-ea)*(fa*ia-ga*ha),ja=(R-V)*(ga-ia)-(Y-ea)*(fa-ha);if(ja){var $=$/ja,aa=aa/ja,ja=+$.toFixed(2),ba=+aa.toFixed(2);R=ja<+X(R,V).toFixed(2)||ja>+W(R,V).toFixed(2)||ja<+X(fa,ha).toFixed(2)||
ja>+W(fa,ha).toFixed(2)||ba<+X(Y,ea).toFixed(2)||ba>+W(Y,ea).toFixed(2)||ba<+X(ga,ia).toFixed(2)||ba>+W(ga,ia).toFixed(2)?void 0:{x:$,y:aa}}else R=void 0}R&&F[R.x.toFixed(4)]!=R.y.toFixed(4)&&(F[R.x.toFixed(4)]=R.y.toFixed(4),Q=Q.t+Z((R[N]-Q[N])/(L[N]-Q[N]))*(L.t-Q.t),B=B.t+Z((R[S]-B[S])/(C[S]-B[S]))*(C.t-B.t),0<=Q&&1>=Q&&0<=B&&1>=B&&(z?M++:M.push({x:R.x,y:R.y,t1:Q,t2:B})))}z=M}else z=z?0:[];if(d)t+=z;else{H=0;for(J=z.length;H<J;H++)z[H].segment1=w,z[H].segment2=G,z[H].bez1=x,z[H].bez2=c;t=t.concat(z)}}}return t}
function r(a){var b=A(a);if(b.bbox)return C(b.bbox);if(!a)return w();a=I(a);for(var d=0,e=0,h=[],f=[],l,n=0,k=a.length;n<k;n++)l=a[n],"M"==l[0]?(d=l[1],e=l[2],h.push(d),f.push(e)):(d=U(d,e,l[1],l[2],l[3],l[4],l[5],l[6]),h=h.concat(d.min.x,d.max.x),f=f.concat(d.min.y,d.max.y),d=l[5],e=l[6]);a=X.apply(0,h);l=X.apply(0,f);h=W.apply(0,h);f=W.apply(0,f);f=w(a,l,h-a,f-l);b.bbox=C(f);return f}function s(a,b,d,e,h){if(h)return[["M",+a+ +h,b],["l",d-2*h,0],["a",h,h,0,0,1,h,h],["l",0,e-2*h],["a",h,h,0,0,1,
-h,h],["l",2*h-d,0],["a",h,h,0,0,1,-h,-h],["l",0,2*h-e],["a",h,h,0,0,1,h,-h],["z"] ];a=[["M",a,b],["l",d,0],["l",0,e],["l",-d,0],["z"] ];a.toString=z;return a}function x(a,b,d,e,h){null==h&&null==e&&(e=d);a=+a;b=+b;d=+d;e=+e;if(null!=h){var f=Math.PI/180,l=a+d*Math.cos(-e*f);a+=d*Math.cos(-h*f);var n=b+d*Math.sin(-e*f);b+=d*Math.sin(-h*f);d=[["M",l,n],["A",d,d,0,+(180<h-e),0,a,b] ]}else d=[["M",a,b],["m",0,-e],["a",d,e,0,1,1,0,2*e],["a",d,e,0,1,1,0,-2*e],["z"] ];d.toString=z;return d}function G(b){var e=
A(b);if(e.abs)return d(e.abs);Q(b,"array")&&Q(b&&b[0],"array")||(b=a.parsePathString(b));if(!b||!b.length)return[["M",0,0] ];var h=[],f=0,l=0,n=0,k=0,p=0;"M"==b[0][0]&&(f=+b[0][1],l=+b[0][2],n=f,k=l,p++,h[0]=["M",f,l]);for(var q=3==b.length&&"M"==b[0][0]&&"R"==b[1][0].toUpperCase()&&"Z"==b[2][0].toUpperCase(),s,r,w=p,c=b.length;w<c;w++){h.push(s=[]);r=b[w];p=r[0];if(p!=p.toUpperCase())switch(s[0]=p.toUpperCase(),s[0]){case "A":s[1]=r[1];s[2]=r[2];s[3]=r[3];s[4]=r[4];s[5]=r[5];s[6]=+r[6]+f;s[7]=+r[7]+
l;break;case "V":s[1]=+r[1]+l;break;case "H":s[1]=+r[1]+f;break;case "R":for(var t=[f,l].concat(r.slice(1)),u=2,v=t.length;u<v;u++)t[u]=+t[u]+f,t[++u]=+t[u]+l;h.pop();h=h.concat(P(t,q));break;case "O":h.pop();t=x(f,l,r[1],r[2]);t.push(t[0]);h=h.concat(t);break;case "U":h.pop();h=h.concat(x(f,l,r[1],r[2],r[3]));s=["U"].concat(h[h.length-1].slice(-2));break;case "M":n=+r[1]+f,k=+r[2]+l;default:for(u=1,v=r.length;u<v;u++)s[u]=+r[u]+(u%2?f:l)}else if("R"==p)t=[f,l].concat(r.slice(1)),h.pop(),h=h.concat(P(t,
q)),s=["R"].concat(r.slice(-2));else if("O"==p)h.pop(),t=x(f,l,r[1],r[2]),t.push(t[0]),h=h.concat(t);else if("U"==p)h.pop(),h=h.concat(x(f,l,r[1],r[2],r[3])),s=["U"].concat(h[h.length-1].slice(-2));else for(t=0,u=r.length;t<u;t++)s[t]=r[t];p=p.toUpperCase();if("O"!=p)switch(s[0]){case "Z":f=+n;l=+k;break;case "H":f=s[1];break;case "V":l=s[1];break;case "M":n=s[s.length-2],k=s[s.length-1];default:f=s[s.length-2],l=s[s.length-1]}}h.toString=z;e.abs=d(h);return h}function h(a,b,d,e){return[a,b,d,e,d,
e]}function J(a,b,d,e,h,f){var l=1/3,n=2/3;return[l*a+n*d,l*b+n*e,l*h+n*d,l*f+n*e,h,f]}function K(b,d,e,h,f,l,n,k,p,s){var r=120*S/180,q=S/180*(+f||0),c=[],t,x=a._.cacher(function(a,b,c){var d=a*F.cos(c)-b*F.sin(c);a=a*F.sin(c)+b*F.cos(c);return{x:d,y:a}});if(s)v=s[0],t=s[1],l=s[2],u=s[3];else{t=x(b,d,-q);b=t.x;d=t.y;t=x(k,p,-q);k=t.x;p=t.y;F.cos(S/180*f);F.sin(S/180*f);t=(b-k)/2;v=(d-p)/2;u=t*t/(e*e)+v*v/(h*h);1<u&&(u=F.sqrt(u),e*=u,h*=u);var u=e*e,w=h*h,u=(l==n?-1:1)*F.sqrt(Z((u*w-u*v*v-w*t*t)/
(u*v*v+w*t*t)));l=u*e*v/h+(b+k)/2;var u=u*-h*t/e+(d+p)/2,v=F.asin(((d-u)/h).toFixed(9));t=F.asin(((p-u)/h).toFixed(9));v=b<l?S-v:v;t=k<l?S-t:t;0>v&&(v=2*S+v);0>t&&(t=2*S+t);n&&v>t&&(v-=2*S);!n&&t>v&&(t-=2*S)}if(Z(t-v)>r){var c=t,w=k,G=p;t=v+r*(n&&t>v?1:-1);k=l+e*F.cos(t);p=u+h*F.sin(t);c=K(k,p,e,h,f,0,n,w,G,[t,c,l,u])}l=t-v;f=F.cos(v);r=F.sin(v);n=F.cos(t);t=F.sin(t);l=F.tan(l/4);e=4/3*e*l;l*=4/3*h;h=[b,d];b=[b+e*r,d-l*f];d=[k+e*t,p-l*n];k=[k,p];b[0]=2*h[0]-b[0];b[1]=2*h[1]-b[1];if(s)return[b,d,k].concat(c);
c=[b,d,k].concat(c).join().split(",");s=[];k=0;for(p=c.length;k<p;k++)s[k]=k%2?x(c[k-1],c[k],q).y:x(c[k],c[k+1],q).x;return s}function U(a,b,d,e,h,f,l,k){for(var n=[],p=[[],[] ],s,r,c,t,q=0;2>q;++q)0==q?(r=6*a-12*d+6*h,s=-3*a+9*d-9*h+3*l,c=3*d-3*a):(r=6*b-12*e+6*f,s=-3*b+9*e-9*f+3*k,c=3*e-3*b),1E-12>Z(s)?1E-12>Z(r)||(s=-c/r,0<s&&1>s&&n.push(s)):(t=r*r-4*c*s,c=F.sqrt(t),0>t||(t=(-r+c)/(2*s),0<t&&1>t&&n.push(t),s=(-r-c)/(2*s),0<s&&1>s&&n.push(s)));for(r=q=n.length;q--;)s=n[q],c=1-s,p[0][q]=c*c*c*a+3*
c*c*s*d+3*c*s*s*h+s*s*s*l,p[1][q]=c*c*c*b+3*c*c*s*e+3*c*s*s*f+s*s*s*k;p[0][r]=a;p[1][r]=b;p[0][r+1]=l;p[1][r+1]=k;p[0].length=p[1].length=r+2;return{min:{x:X.apply(0,p[0]),y:X.apply(0,p[1])},max:{x:W.apply(0,p[0]),y:W.apply(0,p[1])}}}function I(a,b){var e=!b&&A(a);if(!b&&e.curve)return d(e.curve);var f=G(a),l=b&&G(b),n={x:0,y:0,bx:0,by:0,X:0,Y:0,qx:null,qy:null},k={x:0,y:0,bx:0,by:0,X:0,Y:0,qx:null,qy:null},p=function(a,b,c){if(!a)return["C",b.x,b.y,b.x,b.y,b.x,b.y];a[0]in{T:1,Q:1}||(b.qx=b.qy=null);
switch(a[0]){case "M":b.X=a[1];b.Y=a[2];break;case "A":a=["C"].concat(K.apply(0,[b.x,b.y].concat(a.slice(1))));break;case "S":"C"==c||"S"==c?(c=2*b.x-b.bx,b=2*b.y-b.by):(c=b.x,b=b.y);a=["C",c,b].concat(a.slice(1));break;case "T":"Q"==c||"T"==c?(b.qx=2*b.x-b.qx,b.qy=2*b.y-b.qy):(b.qx=b.x,b.qy=b.y);a=["C"].concat(J(b.x,b.y,b.qx,b.qy,a[1],a[2]));break;case "Q":b.qx=a[1];b.qy=a[2];a=["C"].concat(J(b.x,b.y,a[1],a[2],a[3],a[4]));break;case "L":a=["C"].concat(h(b.x,b.y,a[1],a[2]));break;case "H":a=["C"].concat(h(b.x,
b.y,a[1],b.y));break;case "V":a=["C"].concat(h(b.x,b.y,b.x,a[1]));break;case "Z":a=["C"].concat(h(b.x,b.y,b.X,b.Y))}return a},s=function(a,b){if(7<a[b].length){a[b].shift();for(var c=a[b];c.length;)q[b]="A",l&&(u[b]="A"),a.splice(b++,0,["C"].concat(c.splice(0,6)));a.splice(b,1);v=W(f.length,l&&l.length||0)}},r=function(a,b,c,d,e){a&&b&&"M"==a[e][0]&&"M"!=b[e][0]&&(b.splice(e,0,["M",d.x,d.y]),c.bx=0,c.by=0,c.x=a[e][1],c.y=a[e][2],v=W(f.length,l&&l.length||0))},q=[],u=[],c="",t="",x=0,v=W(f.length,
l&&l.length||0);for(;x<v;x++){f[x]&&(c=f[x][0]);"C"!=c&&(q[x]=c,x&&(t=q[x-1]));f[x]=p(f[x],n,t);"A"!=q[x]&&"C"==c&&(q[x]="C");s(f,x);l&&(l[x]&&(c=l[x][0]),"C"!=c&&(u[x]=c,x&&(t=u[x-1])),l[x]=p(l[x],k,t),"A"!=u[x]&&"C"==c&&(u[x]="C"),s(l,x));r(f,l,n,k,x);r(l,f,k,n,x);var w=f[x],z=l&&l[x],y=w.length,U=l&&z.length;n.x=w[y-2];n.y=w[y-1];n.bx=$(w[y-4])||n.x;n.by=$(w[y-3])||n.y;k.bx=l&&($(z[U-4])||k.x);k.by=l&&($(z[U-3])||k.y);k.x=l&&z[U-2];k.y=l&&z[U-1]}l||(e.curve=d(f));return l?[f,l]:f}function P(a,
b){for(var d=[],e=0,h=a.length;h-2*!b>e;e+=2){var f=[{x:+a[e-2],y:+a[e-1]},{x:+a[e],y:+a[e+1]},{x:+a[e+2],y:+a[e+3]},{x:+a[e+4],y:+a[e+5]}];b?e?h-4==e?f[3]={x:+a[0],y:+a[1]}:h-2==e&&(f[2]={x:+a[0],y:+a[1]},f[3]={x:+a[2],y:+a[3]}):f[0]={x:+a[h-2],y:+a[h-1]}:h-4==e?f[3]=f[2]:e||(f[0]={x:+a[e],y:+a[e+1]});d.push(["C",(-f[0].x+6*f[1].x+f[2].x)/6,(-f[0].y+6*f[1].y+f[2].y)/6,(f[1].x+6*f[2].x-f[3].x)/6,(f[1].y+6*f[2].y-f[3].y)/6,f[2].x,f[2].y])}return d}y=k.prototype;var Q=a.is,C=a._.clone,L="hasOwnProperty",
N=/,?([a-z]),?/gi,$=parseFloat,F=Math,S=F.PI,X=F.min,W=F.max,ma=F.pow,Z=F.abs;M=n(1);var na=n(),ba=n(0,1),V=a._unit2px;a.path=A;a.path.getTotalLength=M;a.path.getPointAtLength=na;a.path.getSubpath=function(a,b,d){if(1E-6>this.getTotalLength(a)-d)return ba(a,b).end;a=ba(a,d,1);return b?ba(a,b).end:a};y.getTotalLength=function(){if(this.node.getTotalLength)return this.node.getTotalLength()};y.getPointAtLength=function(a){return na(this.attr("d"),a)};y.getSubpath=function(b,d){return a.path.getSubpath(this.attr("d"),
b,d)};a._.box=w;a.path.findDotsAtSegment=u;a.path.bezierBBox=p;a.path.isPointInsideBBox=b;a.path.isBBoxIntersect=q;a.path.intersection=function(a,b){return l(a,b)};a.path.intersectionNumber=function(a,b){return l(a,b,1)};a.path.isPointInside=function(a,d,e){var h=r(a);return b(h,d,e)&&1==l(a,[["M",d,e],["H",h.x2+10] ],1)%2};a.path.getBBox=r;a.path.get={path:function(a){return a.attr("path")},circle:function(a){a=V(a);return x(a.cx,a.cy,a.r)},ellipse:function(a){a=V(a);return x(a.cx||0,a.cy||0,a.rx,
a.ry)},rect:function(a){a=V(a);return s(a.x||0,a.y||0,a.width,a.height,a.rx,a.ry)},image:function(a){a=V(a);return s(a.x||0,a.y||0,a.width,a.height)},line:function(a){return"M"+[a.attr("x1")||0,a.attr("y1")||0,a.attr("x2"),a.attr("y2")]},polyline:function(a){return"M"+a.attr("points")},polygon:function(a){return"M"+a.attr("points")+"z"},deflt:function(a){a=a.node.getBBox();return s(a.x,a.y,a.width,a.height)}};a.path.toRelative=function(b){var e=A(b),h=String.prototype.toLowerCase;if(e.rel)return d(e.rel);
a.is(b,"array")&&a.is(b&&b[0],"array")||(b=a.parsePathString(b));var f=[],l=0,n=0,k=0,p=0,s=0;"M"==b[0][0]&&(l=b[0][1],n=b[0][2],k=l,p=n,s++,f.push(["M",l,n]));for(var r=b.length;s<r;s++){var q=f[s]=[],x=b[s];if(x[0]!=h.call(x[0]))switch(q[0]=h.call(x[0]),q[0]){case "a":q[1]=x[1];q[2]=x[2];q[3]=x[3];q[4]=x[4];q[5]=x[5];q[6]=+(x[6]-l).toFixed(3);q[7]=+(x[7]-n).toFixed(3);break;case "v":q[1]=+(x[1]-n).toFixed(3);break;case "m":k=x[1],p=x[2];default:for(var c=1,t=x.length;c<t;c++)q[c]=+(x[c]-(c%2?l:
n)).toFixed(3)}else for(f[s]=[],"m"==x[0]&&(k=x[1]+l,p=x[2]+n),q=0,c=x.length;q<c;q++)f[s][q]=x[q];x=f[s].length;switch(f[s][0]){case "z":l=k;n=p;break;case "h":l+=+f[s][x-1];break;case "v":n+=+f[s][x-1];break;default:l+=+f[s][x-2],n+=+f[s][x-1]}}f.toString=z;e.rel=d(f);return f};a.path.toAbsolute=G;a.path.toCubic=I;a.path.map=function(a,b){if(!b)return a;var d,e,h,f,l,n,k;a=I(a);h=0;for(l=a.length;h<l;h++)for(k=a[h],f=1,n=k.length;f<n;f+=2)d=b.x(k[f],k[f+1]),e=b.y(k[f],k[f+1]),k[f]=d,k[f+1]=e;return a};
a.path.toString=z;a.path.clone=d});C.plugin(function(a,v,y,C){var A=Math.max,w=Math.min,z=function(a){this.items=[];this.bindings={};this.length=0;this.type="set";if(a)for(var f=0,n=a.length;f<n;f++)a[f]&&(this[this.items.length]=this.items[this.items.length]=a[f],this.length++)};v=z.prototype;v.push=function(){for(var a,f,n=0,k=arguments.length;n<k;n++)if(a=arguments[n])f=this.items.length,this[f]=this.items[f]=a,this.length++;return this};v.pop=function(){this.length&&delete this[this.length--];
return this.items.pop()};v.forEach=function(a,f){for(var n=0,k=this.items.length;n<k&&!1!==a.call(f,this.items[n],n);n++);return this};v.animate=function(d,f,n,u){"function"!=typeof n||n.length||(u=n,n=L.linear);d instanceof a._.Animation&&(u=d.callback,n=d.easing,f=n.dur,d=d.attr);var p=arguments;if(a.is(d,"array")&&a.is(p[p.length-1],"array"))var b=!0;var q,e=function(){q?this.b=q:q=this.b},l=0,r=u&&function(){l++==this.length&&u.call(this)};return this.forEach(function(a,l){k.once("snap.animcreated."+
a.id,e);b?p[l]&&a.animate.apply(a,p[l]):a.animate(d,f,n,r)})};v.remove=function(){for(;this.length;)this.pop().remove();return this};v.bind=function(a,f,k){var u={};if("function"==typeof f)this.bindings[a]=f;else{var p=k||a;this.bindings[a]=function(a){u[p]=a;f.attr(u)}}return this};v.attr=function(a){var f={},k;for(k in a)if(this.bindings[k])this.bindings[k](a[k]);else f[k]=a[k];a=0;for(k=this.items.length;a<k;a++)this.items[a].attr(f);return this};v.clear=function(){for(;this.length;)this.pop()};
v.splice=function(a,f,k){a=0>a?A(this.length+a,0):a;f=A(0,w(this.length-a,f));var u=[],p=[],b=[],q;for(q=2;q<arguments.length;q++)b.push(arguments[q]);for(q=0;q<f;q++)p.push(this[a+q]);for(;q<this.length-a;q++)u.push(this[a+q]);var e=b.length;for(q=0;q<e+u.length;q++)this.items[a+q]=this[a+q]=q<e?b[q]:u[q-e];for(q=this.items.length=this.length-=f-e;this[q];)delete this[q++];return new z(p)};v.exclude=function(a){for(var f=0,k=this.length;f<k;f++)if(this[f]==a)return this.splice(f,1),!0;return!1};
v.insertAfter=function(a){for(var f=this.items.length;f--;)this.items[f].insertAfter(a);return this};v.getBBox=function(){for(var a=[],f=[],k=[],u=[],p=this.items.length;p--;)if(!this.items[p].removed){var b=this.items[p].getBBox();a.push(b.x);f.push(b.y);k.push(b.x+b.width);u.push(b.y+b.height)}a=w.apply(0,a);f=w.apply(0,f);k=A.apply(0,k);u=A.apply(0,u);return{x:a,y:f,x2:k,y2:u,width:k-a,height:u-f,cx:a+(k-a)/2,cy:f+(u-f)/2}};v.clone=function(a){a=new z;for(var f=0,k=this.items.length;f<k;f++)a.push(this.items[f].clone());
return a};v.toString=function(){return"Snap\u2018s set"};v.type="set";a.set=function(){var a=new z;arguments.length&&a.push.apply(a,Array.prototype.slice.call(arguments,0));return a}});C.plugin(function(a,v,y,C){function A(a){var b=a[0];switch(b.toLowerCase()){case "t":return[b,0,0];case "m":return[b,1,0,0,1,0,0];case "r":return 4==a.length?[b,0,a[2],a[3] ]:[b,0];case "s":return 5==a.length?[b,1,1,a[3],a[4] ]:3==a.length?[b,1,1]:[b,1]}}function w(b,d,f){d=q(d).replace(/\.{3}|\u2026/g,b);b=a.parseTransformString(b)||
[];d=a.parseTransformString(d)||[];for(var k=Math.max(b.length,d.length),p=[],v=[],h=0,w,z,y,I;h<k;h++){y=b[h]||A(d[h]);I=d[h]||A(y);if(y[0]!=I[0]||"r"==y[0].toLowerCase()&&(y[2]!=I[2]||y[3]!=I[3])||"s"==y[0].toLowerCase()&&(y[3]!=I[3]||y[4]!=I[4])){b=a._.transform2matrix(b,f());d=a._.transform2matrix(d,f());p=[["m",b.a,b.b,b.c,b.d,b.e,b.f] ];v=[["m",d.a,d.b,d.c,d.d,d.e,d.f] ];break}p[h]=[];v[h]=[];w=0;for(z=Math.max(y.length,I.length);w<z;w++)w in y&&(p[h][w]=y[w]),w in I&&(v[h][w]=I[w])}return{from:u(p),
to:u(v),f:n(p)}}function z(a){return a}function d(a){return function(b){return+b.toFixed(3)+a}}function f(b){return a.rgb(b[0],b[1],b[2])}function n(a){var b=0,d,f,k,n,h,p,q=[];d=0;for(f=a.length;d<f;d++){h="[";p=['"'+a[d][0]+'"'];k=1;for(n=a[d].length;k<n;k++)p[k]="val["+b++ +"]";h+=p+"]";q[d]=h}return Function("val","return Snap.path.toString.call(["+q+"])")}function u(a){for(var b=[],d=0,f=a.length;d<f;d++)for(var k=1,n=a[d].length;k<n;k++)b.push(a[d][k]);return b}var p={},b=/[a-z]+$/i,q=String;
p.stroke=p.fill="colour";v.prototype.equal=function(a,b){return k("snap.util.equal",this,a,b).firstDefined()};k.on("snap.util.equal",function(e,k){var r,s;r=q(this.attr(e)||"");var x=this;if(r==+r&&k==+k)return{from:+r,to:+k,f:z};if("colour"==p[e])return r=a.color(r),s=a.color(k),{from:[r.r,r.g,r.b,r.opacity],to:[s.r,s.g,s.b,s.opacity],f:f};if("transform"==e||"gradientTransform"==e||"patternTransform"==e)return k instanceof a.Matrix&&(k=k.toTransformString()),a._.rgTransform.test(k)||(k=a._.svgTransform2string(k)),
w(r,k,function(){return x.getBBox(1)});if("d"==e||"path"==e)return r=a.path.toCubic(r,k),{from:u(r[0]),to:u(r[1]),f:n(r[0])};if("points"==e)return r=q(r).split(a._.separator),s=q(k).split(a._.separator),{from:r,to:s,f:function(a){return a}};aUnit=r.match(b);s=q(k).match(b);return aUnit&&aUnit==s?{from:parseFloat(r),to:parseFloat(k),f:d(aUnit)}:{from:this.asPX(e),to:this.asPX(e,k),f:z}})});C.plugin(function(a,v,y,C){var A=v.prototype,w="createTouch"in C.doc;v="click dblclick mousedown mousemove mouseout mouseover mouseup touchstart touchmove touchend touchcancel".split(" ");
var z={mousedown:"touchstart",mousemove:"touchmove",mouseup:"touchend"},d=function(a,b){var d="y"==a?"scrollTop":"scrollLeft",e=b&&b.node?b.node.ownerDocument:C.doc;return e[d in e.documentElement?"documentElement":"body"][d]},f=function(){this.returnValue=!1},n=function(){return this.originalEvent.preventDefault()},u=function(){this.cancelBubble=!0},p=function(){return this.originalEvent.stopPropagation()},b=function(){if(C.doc.addEventListener)return function(a,b,e,f){var k=w&&z[b]?z[b]:b,l=function(k){var l=
d("y",f),q=d("x",f);if(w&&z.hasOwnProperty(b))for(var r=0,u=k.targetTouches&&k.targetTouches.length;r<u;r++)if(k.targetTouches[r].target==a||a.contains(k.targetTouches[r].target)){u=k;k=k.targetTouches[r];k.originalEvent=u;k.preventDefault=n;k.stopPropagation=p;break}return e.call(f,k,k.clientX+q,k.clientY+l)};b!==k&&a.addEventListener(b,l,!1);a.addEventListener(k,l,!1);return function(){b!==k&&a.removeEventListener(b,l,!1);a.removeEventListener(k,l,!1);return!0}};if(C.doc.attachEvent)return function(a,
b,e,h){var k=function(a){a=a||h.node.ownerDocument.window.event;var b=d("y",h),k=d("x",h),k=a.clientX+k,b=a.clientY+b;a.preventDefault=a.preventDefault||f;a.stopPropagation=a.stopPropagation||u;return e.call(h,a,k,b)};a.attachEvent("on"+b,k);return function(){a.detachEvent("on"+b,k);return!0}}}(),q=[],e=function(a){for(var b=a.clientX,e=a.clientY,f=d("y"),l=d("x"),n,p=q.length;p--;){n=q[p];if(w)for(var r=a.touches&&a.touches.length,u;r--;){if(u=a.touches[r],u.identifier==n.el._drag.id||n.el.node.contains(u.target)){b=
u.clientX;e=u.clientY;(a.originalEvent?a.originalEvent:a).preventDefault();break}}else a.preventDefault();b+=l;e+=f;k("snap.drag.move."+n.el.id,n.move_scope||n.el,b-n.el._drag.x,e-n.el._drag.y,b,e,a)}},l=function(b){a.unmousemove(e).unmouseup(l);for(var d=q.length,f;d--;)f=q[d],f.el._drag={},k("snap.drag.end."+f.el.id,f.end_scope||f.start_scope||f.move_scope||f.el,b);q=[]};for(y=v.length;y--;)(function(d){a[d]=A[d]=function(e,f){a.is(e,"function")&&(this.events=this.events||[],this.events.push({name:d,
f:e,unbind:b(this.node||document,d,e,f||this)}));return this};a["un"+d]=A["un"+d]=function(a){for(var b=this.events||[],e=b.length;e--;)if(b[e].name==d&&(b[e].f==a||!a)){b[e].unbind();b.splice(e,1);!b.length&&delete this.events;break}return this}})(v[y]);A.hover=function(a,b,d,e){return this.mouseover(a,d).mouseout(b,e||d)};A.unhover=function(a,b){return this.unmouseover(a).unmouseout(b)};var r=[];A.drag=function(b,d,f,h,n,p){function u(r,v,w){(r.originalEvent||r).preventDefault();this._drag.x=v;
this._drag.y=w;this._drag.id=r.identifier;!q.length&&a.mousemove(e).mouseup(l);q.push({el:this,move_scope:h,start_scope:n,end_scope:p});d&&k.on("snap.drag.start."+this.id,d);b&&k.on("snap.drag.move."+this.id,b);f&&k.on("snap.drag.end."+this.id,f);k("snap.drag.start."+this.id,n||h||this,v,w,r)}if(!arguments.length){var v;return this.drag(function(a,b){this.attr({transform:v+(v?"T":"t")+[a,b]})},function(){v=this.transform().local})}this._drag={};r.push({el:this,start:u});this.mousedown(u);return this};
A.undrag=function(){for(var b=r.length;b--;)r[b].el==this&&(this.unmousedown(r[b].start),r.splice(b,1),k.unbind("snap.drag.*."+this.id));!r.length&&a.unmousemove(e).unmouseup(l);return this}});C.plugin(function(a,v,y,C){y=y.prototype;var A=/^\s*url\((.+)\)/,w=String,z=a._.$;a.filter={};y.filter=function(d){var f=this;"svg"!=f.type&&(f=f.paper);d=a.parse(w(d));var k=a._.id(),u=z("filter");z(u,{id:k,filterUnits:"userSpaceOnUse"});u.appendChild(d.node);f.defs.appendChild(u);return new v(u)};k.on("snap.util.getattr.filter",
function(){k.stop();var d=z(this.node,"filter");if(d)return(d=w(d).match(A))&&a.select(d[1])});k.on("snap.util.attr.filter",function(d){if(d instanceof v&&"filter"==d.type){k.stop();var f=d.node.id;f||(z(d.node,{id:d.id}),f=d.id);z(this.node,{filter:a.url(f)})}d&&"none"!=d||(k.stop(),this.node.removeAttribute("filter"))});a.filter.blur=function(d,f){null==d&&(d=2);return a.format('<feGaussianBlur stdDeviation="{def}"/>',{def:null==f?d:[d,f]})};a.filter.blur.toString=function(){return this()};a.filter.shadow=
function(d,f,k,u,p){"string"==typeof k&&(p=u=k,k=4);"string"!=typeof u&&(p=u,u="#000");null==k&&(k=4);null==p&&(p=1);null==d&&(d=0,f=2);null==f&&(f=d);u=a.color(u||"#000");return a.format('<feGaussianBlur in="SourceAlpha" stdDeviation="{blur}"/><feOffset dx="{dx}" dy="{dy}" result="offsetblur"/><feFlood flood-color="{color}"/><feComposite in2="offsetblur" operator="in"/><feComponentTransfer><feFuncA type="linear" slope="{opacity}"/></feComponentTransfer><feMerge><feMergeNode/><feMergeNode in="SourceGraphic"/></feMerge>',
{color:u,dx:d,dy:f,blur:k,opacity:p})};a.filter.shadow.toString=function(){return this()};a.filter.grayscale=function(d){null==d&&(d=1);return a.format('<feColorMatrix type="matrix" values="{a} {b} {c} 0 0 {d} {e} {f} 0 0 {g} {b} {h} 0 0 0 0 0 1 0"/>',{a:0.2126+0.7874*(1-d),b:0.7152-0.7152*(1-d),c:0.0722-0.0722*(1-d),d:0.2126-0.2126*(1-d),e:0.7152+0.2848*(1-d),f:0.0722-0.0722*(1-d),g:0.2126-0.2126*(1-d),h:0.0722+0.9278*(1-d)})};a.filter.grayscale.toString=function(){return this()};a.filter.sepia=
function(d){null==d&&(d=1);return a.format('<feColorMatrix type="matrix" values="{a} {b} {c} 0 0 {d} {e} {f} 0 0 {g} {h} {i} 0 0 0 0 0 1 0"/>',{a:0.393+0.607*(1-d),b:0.769-0.769*(1-d),c:0.189-0.189*(1-d),d:0.349-0.349*(1-d),e:0.686+0.314*(1-d),f:0.168-0.168*(1-d),g:0.272-0.272*(1-d),h:0.534-0.534*(1-d),i:0.131+0.869*(1-d)})};a.filter.sepia.toString=function(){return this()};a.filter.saturate=function(d){null==d&&(d=1);return a.format('<feColorMatrix type="saturate" values="{amount}"/>',{amount:1-
d})};a.filter.saturate.toString=function(){return this()};a.filter.hueRotate=function(d){return a.format('<feColorMatrix type="hueRotate" values="{angle}"/>',{angle:d||0})};a.filter.hueRotate.toString=function(){return this()};a.filter.invert=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="table" tableValues="{amount} {amount2}"/><feFuncG type="table" tableValues="{amount} {amount2}"/><feFuncB type="table" tableValues="{amount} {amount2}"/></feComponentTransfer>',{amount:d,
amount2:1-d})};a.filter.invert.toString=function(){return this()};a.filter.brightness=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="linear" slope="{amount}"/><feFuncG type="linear" slope="{amount}"/><feFuncB type="linear" slope="{amount}"/></feComponentTransfer>',{amount:d})};a.filter.brightness.toString=function(){return this()};a.filter.contrast=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="linear" slope="{amount}" intercept="{amount2}"/><feFuncG type="linear" slope="{amount}" intercept="{amount2}"/><feFuncB type="linear" slope="{amount}" intercept="{amount2}"/></feComponentTransfer>',
{amount:d,amount2:0.5-d/2})};a.filter.contrast.toString=function(){return this()}});return C});

]]> </script>
</svg>


<pre class='language-julia'><code class='language-julia'>mat(control(3, 3, 1=&gt;shift(π/4)))</code></pre>
<pre id='var-hash271816' class='code-output documenter-example-output'>8×8 Diagonal{ComplexF64, Vector{ComplexF64}}:
 1.0+0.0im      ⋅          ⋅          ⋅      …      ⋅               ⋅    
     ⋅      1.0+0.0im      ⋅          ⋅             ⋅               ⋅    
     ⋅          ⋅      1.0+0.0im      ⋅             ⋅               ⋅    
     ⋅          ⋅          ⋅      1.0+0.0im         ⋅               ⋅    
     ⋅          ⋅          ⋅          ⋅             ⋅               ⋅    
     ⋅          ⋅          ⋅          ⋅      …      ⋅               ⋅    
     ⋅          ⋅          ⋅          ⋅         1.0+0.0im           ⋅    
     ⋅          ⋅          ⋅          ⋅             ⋅      0.707107+0.707107im</pre>


<div class="markdown"><p><strong>Yao blocks are fast</strong></p>
</div>

<pre class='language-julia'><code class='language-julia'>let
    x = @bind run_benchmark_mat CheckBox()
    md" $x run the benchmark"
end</code></pre>
<div class="markdown"><p><bond def="run_benchmark_mat" unique_id="zyO2bVSxjc0Q"><input type="checkbox"></bond> run the benchmark</p>
</div>

<pre class='language-julia'><code class='language-julia'>if run_benchmark_mat
    @benchmark mat(control(16, 7, 3=&gt;shift(π/4)))
end</code></pre>


<pre class='language-julia'><code class='language-julia'>if run_benchmark_mat
    @benchmark kron(
        kron(
            kron(Diagonal(ones(ComplexF64, 2^9)), Diagonal(mat(ConstGate.P1))),
                Diagonal(ones(ComplexF64, 2^3))
        ), 
        kron(mat(ComplexF64, shift(π/4)), Diagonal(ones(ComplexF64, 2^2)))
    ) +
    kron(
        kron(Diagonal(ones(ComplexF64, 2^9)), Diagonal(mat(ConstGate.P0))),
            Diagonal(ones(ComplexF64, 2^6))
    )
end</code></pre>



<div class="markdown"><p><strong>Yao block is general purposed</strong></p>
</div>

<pre class='language-julia'><code class='language-julia'>multi_control = control(10, (2, -6, 4), (9, 1)=&gt;chain(put(2, 1=&gt;X), rot(kron(X, X), 0.3)))</code></pre>
<pre id='var-multi_control' class='code-output documenter-example-output'>nqubits: 10
control(2, ¬6, 4)
└─ (9, 1) chain
   ├─ put on (1)
   │  └─ X
   └─ rot(nqubits: 2
kron
├─ 1=>X
└─ 2=>X, 0.3)
</pre>

<pre class='language-julia'><code class='language-julia'>vizcircuit(multi_control)</code></pre>
<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:xlink="http://www.w3.org/1999/xlink"
     version="1.2"
     width="85mm" height="187.5mm" viewBox="0 0 85 187.5"
     stroke="none"
     fill="#000000"
     stroke-width="0.3"
     font-size="3.88"

     id="img-34b61313">
<defs>
  <marker id="arrow" markerWidth="15" markerHeight="7" refX="5" refY="3.5" orient="auto" markerUnits="strokeWidth">
    <path d="M0,0 L15,3.5 L0,7 z" stroke="context-stroke" fill="context-stroke"/>
  </marker>
</defs>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-34b61313-1">
  <g transform="translate(55,81.25)">
    <path d="M-16.02,-80 L16.02,-80 16.02,80 -16.02,80  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-34b61313-2">
  <g transform="translate(55,100)">
    <circle cx="0" cy="0" r="0.83" class="primitive"/>
  </g>
</g>
<g stroke-width="0" fill="#000000" id="img-34b61313-3">
  <g transform="translate(55,62.5)">
    <circle cx="0" cy="0" r="0.83" class="primitive"/>
  </g>
</g>
<g stroke-width="0" fill="#000000" id="img-34b61313-4">
  <g transform="translate(55,25)">
    <circle cx="0" cy="0" r="0.83" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-5">
  <g transform="translate(23.12,156.25)">
    <path d="M0,-3.5 L0,3.5  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-6">
  <g transform="translate(23.12,156.25)">
    <path d="M-3.5,0 L3.5,0  z" class="primitive"/>
  </g>
</g>
<g fill="#000000" fill-opacity="0.000" stroke-width="0.35" stroke="#000000" id="img-34b61313-7">
  <g transform="translate(23.12,156.25)">
    <circle cx="0" cy="0" r="3.5" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-34b61313-8">
  <g transform="translate(23.12,100)">
    <circle cx="0" cy="0" r="0.83" class="primitive"/>
  </g>
</g>
<g stroke-width="0" fill="#000000" id="img-34b61313-9">
  <g transform="translate(23.12,62.5)">
    <circle cx="0" cy="0" r="0.83" class="primitive"/>
  </g>
</g>
<g stroke-width="0" fill="#000000" id="img-34b61313-10">
  <g transform="translate(23.12,25)">
    <circle cx="0" cy="0" r="0.83" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-11">
  <g transform="translate(44.38,175)">
    <path fill="none" d="M31.88,0 L-31.88,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-12">
  <g transform="translate(76.25,156.25)">
    <path fill="none" d="M0,0 L0,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-13">
  <g transform="translate(76.25,137.5)">
    <path fill="none" d="M0,0 L0,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-14">
  <g transform="translate(76.25,118.75)">
    <path fill="none" d="M0,0 L0,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-15">
  <g transform="translate(76.25,100)">
    <path fill="none" d="M0,0 L0,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-16">
  <g transform="translate(76.25,81.25)">
    <path fill="none" d="M0,0 L0,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-17">
  <g transform="translate(76.25,62.5)">
    <path fill="none" d="M0,0 L0,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-18">
  <g transform="translate(76.25,43.75)">
    <path fill="none" d="M0,0 L0,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-19">
  <g transform="translate(76.25,25)">
    <path fill="none" d="M0,0 L0,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-20">
  <g transform="translate(76.25,6.25)">
    <path fill="none" d="M0,0 L0,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-21">
  <g transform="translate(55,137.5)">
    <path fill="none" d="M-21.25,0 L21.25,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-22">
  <g transform="translate(55,118.75)">
    <path fill="none" d="M-21.25,0 L21.25,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-23">
  <g transform="translate(55,81.25)">
    <path fill="none" d="M-21.25,0 L21.25,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-24">
  <g transform="translate(55,43.75)">
    <path fill="none" d="M-21.25,0 L21.25,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-25">
  <g transform="translate(73.64,6.25)">
    <path fill="none" d="M-2.61,0 L2.61,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-26">
  <g transform="translate(25.74,6.25)">
    <path fill="none" d="M-13.24,0 L13.24,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-27">
  <g transform="translate(73.64,156.25)">
    <path fill="none" d="M-2.61,0 L2.61,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-28">
  <g transform="translate(36.36,156.25)">
    <path fill="none" d="M-2.61,0 L2.61,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-29">
  <g transform="translate(66.04,100)">
    <path fill="none" d="M-10.21,0 L10.21,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-30">
  <g transform="translate(43.96,100)">
    <path fill="none" d="M-10.21,0 L10.21,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-31">
  <g transform="translate(66.04,62.5)">
    <path fill="none" d="M-10.21,0 L10.21,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-32">
  <g transform="translate(43.96,62.5)">
    <path fill="none" d="M-10.21,0 L10.21,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-33">
  <g transform="translate(66.04,25)">
    <path fill="none" d="M-10.21,0 L10.21,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-34">
  <g transform="translate(43.96,25)">
    <path fill="none" d="M-10.21,0 L10.21,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-35">
  <g transform="translate(55,51.04)">
    <path fill="none" d="M0,-49.79 L0,49.79 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-36">
  <g transform="translate(55,81.25)">
    <path fill="none" d="M0,17.92 L0,-17.92 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-37">
  <g transform="translate(55,43.75)">
    <path fill="none" d="M0,17.92 L0,-17.92 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-38">
  <g transform="translate(23.12,137.5)">
    <path fill="none" d="M-10.62,0 L10.62,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-39">
  <g transform="translate(23.12,118.75)">
    <path fill="none" d="M-10.62,0 L10.62,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-40">
  <g transform="translate(23.12,81.25)">
    <path fill="none" d="M-10.62,0 L10.62,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-41">
  <g transform="translate(23.12,43.75)">
    <path fill="none" d="M-10.62,0 L10.62,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-42">
  <g transform="translate(30.19,156.25)">
    <path fill="none" d="M-3.56,0 L3.56,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-43">
  <g transform="translate(16.06,156.25)">
    <path fill="none" d="M-3.56,0 L3.56,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-44">
  <g transform="translate(28.85,100)">
    <path fill="none" d="M-4.9,0 L4.9,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-45">
  <g transform="translate(17.4,100)">
    <path fill="none" d="M-4.9,0 L4.9,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-46">
  <g transform="translate(28.85,62.5)">
    <path fill="none" d="M-4.9,0 L4.9,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-47">
  <g transform="translate(17.4,62.5)">
    <path fill="none" d="M-4.9,0 L4.9,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-48">
  <g transform="translate(28.85,25)">
    <path fill="none" d="M-4.9,0 L4.9,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-49">
  <g transform="translate(17.4,25)">
    <path fill="none" d="M-4.9,0 L4.9,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-50">
  <g transform="translate(23.12,126.79)">
    <path fill="none" d="M0,25.96 L0,-25.96 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-51">
  <g transform="translate(23.12,81.25)">
    <path fill="none" d="M0,17.92 L0,-17.92 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-34b61313-52">
  <g transform="translate(23.12,43.75)">
    <path fill="none" d="M0,17.92 L0,-17.92 " class="primitive"/>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="3.53" id="img-34b61313-53">
  <g transform="translate(55,81.25)">
    <g class="primitive">
      <text text-anchor="middle" dy="-1.45em">rot(nqubits: 2<tspan x="0" dy="1.1000em">kron</tspan><tspan x="0" dy="1.1000em">├─ 1=>X</tspan><tspan x="0" dy="1.1000em">└─ 2=>X, 0.3)</tspan></text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="5.64" id="img-34b61313-54">
  <g transform="translate(55,100)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em"></text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="5.64" id="img-34b61313-55">
  <g transform="translate(55,62.5)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em"></text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="5.64" id="img-34b61313-56">
  <g transform="translate(55,25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em"></text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="5.64" id="img-34b61313-57">
  <g transform="translate(23.12,156.25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em"></text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="5.64" id="img-34b61313-58">
  <g transform="translate(23.12,100)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em"></text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="5.64" id="img-34b61313-59">
  <g transform="translate(23.12,62.5)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em"></text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="5.64" id="img-34b61313-60">
  <g transform="translate(23.12,25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em"></text>
    </g>
  </g>
</g>
<script> <![CDATA[
(function(N){var k=/[\.\/]/,L=/\s*,\s*/,C=function(a,d){return a-d},a,v,y={n:{}},M=function(){for(var a=0,d=this.length;a<d;a++)if("undefined"!=typeof this[a])return this[a]},A=function(){for(var a=this.length;--a;)if("undefined"!=typeof this[a])return this[a]},w=function(k,d){k=String(k);var f=v,n=Array.prototype.slice.call(arguments,2),u=w.listeners(k),p=0,b,q=[],e={},l=[],r=a;l.firstDefined=M;l.lastDefined=A;a=k;for(var s=v=0,x=u.length;s<x;s++)"zIndex"in u[s]&&(q.push(u[s].zIndex),0>u[s].zIndex&&
(e[u[s].zIndex]=u[s]));for(q.sort(C);0>q[p];)if(b=e[q[p++] ],l.push(b.apply(d,n)),v)return v=f,l;for(s=0;s<x;s++)if(b=u[s],"zIndex"in b)if(b.zIndex==q[p]){l.push(b.apply(d,n));if(v)break;do if(p++,(b=e[q[p] ])&&l.push(b.apply(d,n)),v)break;while(b)}else e[b.zIndex]=b;else if(l.push(b.apply(d,n)),v)break;v=f;a=r;return l};w._events=y;w.listeners=function(a){a=a.split(k);var d=y,f,n,u,p,b,q,e,l=[d],r=[];u=0;for(p=a.length;u<p;u++){e=[];b=0;for(q=l.length;b<q;b++)for(d=l[b].n,f=[d[a[u] ],d["*"] ],n=2;n--;)if(d=
f[n])e.push(d),r=r.concat(d.f||[]);l=e}return r};w.on=function(a,d){a=String(a);if("function"!=typeof d)return function(){};for(var f=a.split(L),n=0,u=f.length;n<u;n++)(function(a){a=a.split(k);for(var b=y,f,e=0,l=a.length;e<l;e++)b=b.n,b=b.hasOwnProperty(a[e])&&b[a[e] ]||(b[a[e] ]={n:{}});b.f=b.f||[];e=0;for(l=b.f.length;e<l;e++)if(b.f[e]==d){f=!0;break}!f&&b.f.push(d)})(f[n]);return function(a){+a==+a&&(d.zIndex=+a)}};w.f=function(a){var d=[].slice.call(arguments,1);return function(){w.apply(null,
[a,null].concat(d).concat([].slice.call(arguments,0)))}};w.stop=function(){v=1};w.nt=function(k){return k?(new RegExp("(?:\\.|\\/|^)"+k+"(?:\\.|\\/|$)")).test(a):a};w.nts=function(){return a.split(k)};w.off=w.unbind=function(a,d){if(a){var f=a.split(L);if(1<f.length)for(var n=0,u=f.length;n<u;n++)w.off(f[n],d);else{for(var f=a.split(k),p,b,q,e,l=[y],n=0,u=f.length;n<u;n++)for(e=0;e<l.length;e+=q.length-2){q=[e,1];p=l[e].n;if("*"!=f[n])p[f[n] ]&&q.push(p[f[n] ]);else for(b in p)p.hasOwnProperty(b)&&
q.push(p[b]);l.splice.apply(l,q)}n=0;for(u=l.length;n<u;n++)for(p=l[n];p.n;){if(d){if(p.f){e=0;for(f=p.f.length;e<f;e++)if(p.f[e]==d){p.f.splice(e,1);break}!p.f.length&&delete p.f}for(b in p.n)if(p.n.hasOwnProperty(b)&&p.n[b].f){q=p.n[b].f;e=0;for(f=q.length;e<f;e++)if(q[e]==d){q.splice(e,1);break}!q.length&&delete p.n[b].f}}else for(b in delete p.f,p.n)p.n.hasOwnProperty(b)&&p.n[b].f&&delete p.n[b].f;p=p.n}}}else w._events=y={n:{}}};w.once=function(a,d){var f=function(){w.unbind(a,f);return d.apply(this,
arguments)};return w.on(a,f)};w.version="0.4.2";w.toString=function(){return"You are running Eve 0.4.2"};"undefined"!=typeof module&&module.exports?module.exports=w:"function"===typeof define&&define.amd?define("eve",[],function(){return w}):N.eve=w})(this);
(function(N,k){"function"===typeof define&&define.amd?define("Snap.svg",["eve"],function(L){return k(N,L)}):k(N,N.eve)})(this,function(N,k){var L=function(a){var k={},y=N.requestAnimationFrame||N.webkitRequestAnimationFrame||N.mozRequestAnimationFrame||N.oRequestAnimationFrame||N.msRequestAnimationFrame||function(a){setTimeout(a,16)},M=Array.isArray||function(a){return a instanceof Array||"[object Array]"==Object.prototype.toString.call(a)},A=0,w="M"+(+new Date).toString(36),z=function(a){if(null==
a)return this.s;var b=this.s-a;this.b+=this.dur*b;this.B+=this.dur*b;this.s=a},d=function(a){if(null==a)return this.spd;this.spd=a},f=function(a){if(null==a)return this.dur;this.s=this.s*a/this.dur;this.dur=a},n=function(){delete k[this.id];this.update();a("mina.stop."+this.id,this)},u=function(){this.pdif||(delete k[this.id],this.update(),this.pdif=this.get()-this.b)},p=function(){this.pdif&&(this.b=this.get()-this.pdif,delete this.pdif,k[this.id]=this)},b=function(){var a;if(M(this.start)){a=[];
for(var b=0,e=this.start.length;b<e;b++)a[b]=+this.start[b]+(this.end[b]-this.start[b])*this.easing(this.s)}else a=+this.start+(this.end-this.start)*this.easing(this.s);this.set(a)},q=function(){var l=0,b;for(b in k)if(k.hasOwnProperty(b)){var e=k[b],f=e.get();l++;e.s=(f-e.b)/(e.dur/e.spd);1<=e.s&&(delete k[b],e.s=1,l--,function(b){setTimeout(function(){a("mina.finish."+b.id,b)})}(e));e.update()}l&&y(q)},e=function(a,r,s,x,G,h,J){a={id:w+(A++).toString(36),start:a,end:r,b:s,s:0,dur:x-s,spd:1,get:G,
set:h,easing:J||e.linear,status:z,speed:d,duration:f,stop:n,pause:u,resume:p,update:b};k[a.id]=a;r=0;for(var K in k)if(k.hasOwnProperty(K)&&(r++,2==r))break;1==r&&y(q);return a};e.time=Date.now||function(){return+new Date};e.getById=function(a){return k[a]||null};e.linear=function(a){return a};e.easeout=function(a){return Math.pow(a,1.7)};e.easein=function(a){return Math.pow(a,0.48)};e.easeinout=function(a){if(1==a)return 1;if(0==a)return 0;var b=0.48-a/1.04,e=Math.sqrt(0.1734+b*b);a=e-b;a=Math.pow(Math.abs(a),
1/3)*(0>a?-1:1);b=-e-b;b=Math.pow(Math.abs(b),1/3)*(0>b?-1:1);a=a+b+0.5;return 3*(1-a)*a*a+a*a*a};e.backin=function(a){return 1==a?1:a*a*(2.70158*a-1.70158)};e.backout=function(a){if(0==a)return 0;a-=1;return a*a*(2.70158*a+1.70158)+1};e.elastic=function(a){return a==!!a?a:Math.pow(2,-10*a)*Math.sin(2*(a-0.075)*Math.PI/0.3)+1};e.bounce=function(a){a<1/2.75?a*=7.5625*a:a<2/2.75?(a-=1.5/2.75,a=7.5625*a*a+0.75):a<2.5/2.75?(a-=2.25/2.75,a=7.5625*a*a+0.9375):(a-=2.625/2.75,a=7.5625*a*a+0.984375);return a};
return N.mina=e}("undefined"==typeof k?function(){}:k),C=function(){function a(c,t){if(c){if(c.tagName)return x(c);if(y(c,"array")&&a.set)return a.set.apply(a,c);if(c instanceof e)return c;if(null==t)return c=G.doc.querySelector(c),x(c)}return new s(null==c?"100%":c,null==t?"100%":t)}function v(c,a){if(a){"#text"==c&&(c=G.doc.createTextNode(a.text||""));"string"==typeof c&&(c=v(c));if("string"==typeof a)return"xlink:"==a.substring(0,6)?c.getAttributeNS(m,a.substring(6)):"xml:"==a.substring(0,4)?c.getAttributeNS(la,
a.substring(4)):c.getAttribute(a);for(var da in a)if(a[h](da)){var b=J(a[da]);b?"xlink:"==da.substring(0,6)?c.setAttributeNS(m,da.substring(6),b):"xml:"==da.substring(0,4)?c.setAttributeNS(la,da.substring(4),b):c.setAttribute(da,b):c.removeAttribute(da)}}else c=G.doc.createElementNS(la,c);return c}function y(c,a){a=J.prototype.toLowerCase.call(a);return"finite"==a?isFinite(c):"array"==a&&(c instanceof Array||Array.isArray&&Array.isArray(c))?!0:"null"==a&&null===c||a==typeof c&&null!==c||"object"==
a&&c===Object(c)||$.call(c).slice(8,-1).toLowerCase()==a}function M(c){if("function"==typeof c||Object(c)!==c)return c;var a=new c.constructor,b;for(b in c)c[h](b)&&(a[b]=M(c[b]));return a}function A(c,a,b){function m(){var e=Array.prototype.slice.call(arguments,0),f=e.join("\u2400"),d=m.cache=m.cache||{},l=m.count=m.count||[];if(d[h](f)){a:for(var e=l,l=f,B=0,H=e.length;B<H;B++)if(e[B]===l){e.push(e.splice(B,1)[0]);break a}return b?b(d[f]):d[f]}1E3<=l.length&&delete d[l.shift()];l.push(f);d[f]=c.apply(a,
e);return b?b(d[f]):d[f]}return m}function w(c,a,b,m,e,f){return null==e?(c-=b,a-=m,c||a?(180*I.atan2(-a,-c)/C+540)%360:0):w(c,a,e,f)-w(b,m,e,f)}function z(c){return c%360*C/180}function d(c){var a=[];c=c.replace(/(?:^|\s)(\w+)\(([^)]+)\)/g,function(c,b,m){m=m.split(/\s*,\s*|\s+/);"rotate"==b&&1==m.length&&m.push(0,0);"scale"==b&&(2<m.length?m=m.slice(0,2):2==m.length&&m.push(0,0),1==m.length&&m.push(m[0],0,0));"skewX"==b?a.push(["m",1,0,I.tan(z(m[0])),1,0,0]):"skewY"==b?a.push(["m",1,I.tan(z(m[0])),
0,1,0,0]):a.push([b.charAt(0)].concat(m));return c});return a}function f(c,t){var b=O(c),m=new a.Matrix;if(b)for(var e=0,f=b.length;e<f;e++){var h=b[e],d=h.length,B=J(h[0]).toLowerCase(),H=h[0]!=B,l=H?m.invert():0,E;"t"==B&&2==d?m.translate(h[1],0):"t"==B&&3==d?H?(d=l.x(0,0),B=l.y(0,0),H=l.x(h[1],h[2]),l=l.y(h[1],h[2]),m.translate(H-d,l-B)):m.translate(h[1],h[2]):"r"==B?2==d?(E=E||t,m.rotate(h[1],E.x+E.width/2,E.y+E.height/2)):4==d&&(H?(H=l.x(h[2],h[3]),l=l.y(h[2],h[3]),m.rotate(h[1],H,l)):m.rotate(h[1],
h[2],h[3])):"s"==B?2==d||3==d?(E=E||t,m.scale(h[1],h[d-1],E.x+E.width/2,E.y+E.height/2)):4==d?H?(H=l.x(h[2],h[3]),l=l.y(h[2],h[3]),m.scale(h[1],h[1],H,l)):m.scale(h[1],h[1],h[2],h[3]):5==d&&(H?(H=l.x(h[3],h[4]),l=l.y(h[3],h[4]),m.scale(h[1],h[2],H,l)):m.scale(h[1],h[2],h[3],h[4])):"m"==B&&7==d&&m.add(h[1],h[2],h[3],h[4],h[5],h[6])}return m}function n(c,t){if(null==t){var m=!0;t="linearGradient"==c.type||"radialGradient"==c.type?c.node.getAttribute("gradientTransform"):"pattern"==c.type?c.node.getAttribute("patternTransform"):
c.node.getAttribute("transform");if(!t)return new a.Matrix;t=d(t)}else t=a._.rgTransform.test(t)?J(t).replace(/\.{3}|\u2026/g,c._.transform||aa):d(t),y(t,"array")&&(t=a.path?a.path.toString.call(t):J(t)),c._.transform=t;var b=f(t,c.getBBox(1));if(m)return b;c.matrix=b}function u(c){c=c.node.ownerSVGElement&&x(c.node.ownerSVGElement)||c.node.parentNode&&x(c.node.parentNode)||a.select("svg")||a(0,0);var t=c.select("defs"),t=null==t?!1:t.node;t||(t=r("defs",c.node).node);return t}function p(c){return c.node.ownerSVGElement&&
x(c.node.ownerSVGElement)||a.select("svg")}function b(c,a,m){function b(c){if(null==c)return aa;if(c==+c)return c;v(B,{width:c});try{return B.getBBox().width}catch(a){return 0}}function h(c){if(null==c)return aa;if(c==+c)return c;v(B,{height:c});try{return B.getBBox().height}catch(a){return 0}}function e(b,B){null==a?d[b]=B(c.attr(b)||0):b==a&&(d=B(null==m?c.attr(b)||0:m))}var f=p(c).node,d={},B=f.querySelector(".svg---mgr");B||(B=v("rect"),v(B,{x:-9E9,y:-9E9,width:10,height:10,"class":"svg---mgr",
fill:"none"}),f.appendChild(B));switch(c.type){case "rect":e("rx",b),e("ry",h);case "image":e("width",b),e("height",h);case "text":e("x",b);e("y",h);break;case "circle":e("cx",b);e("cy",h);e("r",b);break;case "ellipse":e("cx",b);e("cy",h);e("rx",b);e("ry",h);break;case "line":e("x1",b);e("x2",b);e("y1",h);e("y2",h);break;case "marker":e("refX",b);e("markerWidth",b);e("refY",h);e("markerHeight",h);break;case "radialGradient":e("fx",b);e("fy",h);break;case "tspan":e("dx",b);e("dy",h);break;default:e(a,
b)}f.removeChild(B);return d}function q(c){y(c,"array")||(c=Array.prototype.slice.call(arguments,0));for(var a=0,b=0,m=this.node;this[a];)delete this[a++];for(a=0;a<c.length;a++)"set"==c[a].type?c[a].forEach(function(c){m.appendChild(c.node)}):m.appendChild(c[a].node);for(var h=m.childNodes,a=0;a<h.length;a++)this[b++]=x(h[a]);return this}function e(c){if(c.snap in E)return E[c.snap];var a=this.id=V(),b;try{b=c.ownerSVGElement}catch(m){}this.node=c;b&&(this.paper=new s(b));this.type=c.tagName;this.anims=
{};this._={transform:[]};c.snap=a;E[a]=this;"g"==this.type&&(this.add=q);if(this.type in{g:1,mask:1,pattern:1})for(var e in s.prototype)s.prototype[h](e)&&(this[e]=s.prototype[e])}function l(c){this.node=c}function r(c,a){var b=v(c);a.appendChild(b);return x(b)}function s(c,a){var b,m,f,d=s.prototype;if(c&&"svg"==c.tagName){if(c.snap in E)return E[c.snap];var l=c.ownerDocument;b=new e(c);m=c.getElementsByTagName("desc")[0];f=c.getElementsByTagName("defs")[0];m||(m=v("desc"),m.appendChild(l.createTextNode("Created with Snap")),
b.node.appendChild(m));f||(f=v("defs"),b.node.appendChild(f));b.defs=f;for(var ca in d)d[h](ca)&&(b[ca]=d[ca]);b.paper=b.root=b}else b=r("svg",G.doc.body),v(b.node,{height:a,version:1.1,width:c,xmlns:la});return b}function x(c){return!c||c instanceof e||c instanceof l?c:c.tagName&&"svg"==c.tagName.toLowerCase()?new s(c):c.tagName&&"object"==c.tagName.toLowerCase()&&"image/svg+xml"==c.type?new s(c.contentDocument.getElementsByTagName("svg")[0]):new e(c)}a.version="0.3.0";a.toString=function(){return"Snap v"+
this.version};a._={};var G={win:N,doc:N.document};a._.glob=G;var h="hasOwnProperty",J=String,K=parseFloat,U=parseInt,I=Math,P=I.max,Q=I.min,Y=I.abs,C=I.PI,aa="",$=Object.prototype.toString,F=/^\s*((#[a-f\d]{6})|(#[a-f\d]{3})|rgba?\(\s*([\d\.]+%?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+%?(?:\s*,\s*[\d\.]+%?)?)\s*\)|hsba?\(\s*([\d\.]+(?:deg|\xb0|%)?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+(?:%?\s*,\s*[\d\.]+)?%?)\s*\)|hsla?\(\s*([\d\.]+(?:deg|\xb0|%)?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+(?:%?\s*,\s*[\d\.]+)?%?)\s*\))\s*$/i;a._.separator=
RegExp("[,\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]+");var S=RegExp("[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*"),X={hs:1,rg:1},W=RegExp("([a-z])[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029,]*((-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*)+)",
"ig"),ma=RegExp("([rstm])[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029,]*((-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*)+)","ig"),Z=RegExp("(-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?)[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*",
"ig"),na=0,ba="S"+(+new Date).toString(36),V=function(){return ba+(na++).toString(36)},m="http://www.w3.org/1999/xlink",la="http://www.w3.org/2000/svg",E={},ca=a.url=function(c){return"url('#"+c+"')"};a._.$=v;a._.id=V;a.format=function(){var c=/\{([^\}]+)\}/g,a=/(?:(?:^|\.)(.+?)(?=\[|\.|$|\()|\[('|")(.+?)\2\])(\(\))?/g,b=function(c,b,m){var h=m;b.replace(a,function(c,a,b,m,t){a=a||m;h&&(a in h&&(h=h[a]),"function"==typeof h&&t&&(h=h()))});return h=(null==h||h==m?c:h)+""};return function(a,m){return J(a).replace(c,
function(c,a){return b(c,a,m)})}}();a._.clone=M;a._.cacher=A;a.rad=z;a.deg=function(c){return 180*c/C%360};a.angle=w;a.is=y;a.snapTo=function(c,a,b){b=y(b,"finite")?b:10;if(y(c,"array"))for(var m=c.length;m--;){if(Y(c[m]-a)<=b)return c[m]}else{c=+c;m=a%c;if(m<b)return a-m;if(m>c-b)return a-m+c}return a};a.getRGB=A(function(c){if(!c||(c=J(c)).indexOf("-")+1)return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka};if("none"==c)return{r:-1,g:-1,b:-1,hex:"none",toString:ka};!X[h](c.toLowerCase().substring(0,
2))&&"#"!=c.charAt()&&(c=T(c));if(!c)return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka};var b,m,e,f,d;if(c=c.match(F)){c[2]&&(e=U(c[2].substring(5),16),m=U(c[2].substring(3,5),16),b=U(c[2].substring(1,3),16));c[3]&&(e=U((d=c[3].charAt(3))+d,16),m=U((d=c[3].charAt(2))+d,16),b=U((d=c[3].charAt(1))+d,16));c[4]&&(d=c[4].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b*=2.55),m=K(d[1]),"%"==d[1].slice(-1)&&(m*=2.55),e=K(d[2]),"%"==d[2].slice(-1)&&(e*=2.55),"rgba"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),
d[3]&&"%"==d[3].slice(-1)&&(f/=100));if(c[5])return d=c[5].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b/=100),m=K(d[1]),"%"==d[1].slice(-1)&&(m/=100),e=K(d[2]),"%"==d[2].slice(-1)&&(e/=100),"deg"!=d[0].slice(-3)&&"\u00b0"!=d[0].slice(-1)||(b/=360),"hsba"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),d[3]&&"%"==d[3].slice(-1)&&(f/=100),a.hsb2rgb(b,m,e,f);if(c[6])return d=c[6].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b/=100),m=K(d[1]),"%"==d[1].slice(-1)&&(m/=100),e=K(d[2]),"%"==d[2].slice(-1)&&(e/=100),
"deg"!=d[0].slice(-3)&&"\u00b0"!=d[0].slice(-1)||(b/=360),"hsla"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),d[3]&&"%"==d[3].slice(-1)&&(f/=100),a.hsl2rgb(b,m,e,f);b=Q(I.round(b),255);m=Q(I.round(m),255);e=Q(I.round(e),255);f=Q(P(f,0),1);c={r:b,g:m,b:e,toString:ka};c.hex="#"+(16777216|e|m<<8|b<<16).toString(16).slice(1);c.opacity=y(f,"finite")?f:1;return c}return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka}},a);a.hsb=A(function(c,b,m){return a.hsb2rgb(c,b,m).hex});a.hsl=A(function(c,b,m){return a.hsl2rgb(c,
b,m).hex});a.rgb=A(function(c,a,b,m){if(y(m,"finite")){var e=I.round;return"rgba("+[e(c),e(a),e(b),+m.toFixed(2)]+")"}return"#"+(16777216|b|a<<8|c<<16).toString(16).slice(1)});var T=function(c){var a=G.doc.getElementsByTagName("head")[0]||G.doc.getElementsByTagName("svg")[0];T=A(function(c){if("red"==c.toLowerCase())return"rgb(255, 0, 0)";a.style.color="rgb(255, 0, 0)";a.style.color=c;c=G.doc.defaultView.getComputedStyle(a,aa).getPropertyValue("color");return"rgb(255, 0, 0)"==c?null:c});return T(c)},
qa=function(){return"hsb("+[this.h,this.s,this.b]+")"},ra=function(){return"hsl("+[this.h,this.s,this.l]+")"},ka=function(){return 1==this.opacity||null==this.opacity?this.hex:"rgba("+[this.r,this.g,this.b,this.opacity]+")"},D=function(c,b,m){null==b&&y(c,"object")&&"r"in c&&"g"in c&&"b"in c&&(m=c.b,b=c.g,c=c.r);null==b&&y(c,string)&&(m=a.getRGB(c),c=m.r,b=m.g,m=m.b);if(1<c||1<b||1<m)c/=255,b/=255,m/=255;return[c,b,m]},oa=function(c,b,m,e){c=I.round(255*c);b=I.round(255*b);m=I.round(255*m);c={r:c,
g:b,b:m,opacity:y(e,"finite")?e:1,hex:a.rgb(c,b,m),toString:ka};y(e,"finite")&&(c.opacity=e);return c};a.color=function(c){var b;y(c,"object")&&"h"in c&&"s"in c&&"b"in c?(b=a.hsb2rgb(c),c.r=b.r,c.g=b.g,c.b=b.b,c.opacity=1,c.hex=b.hex):y(c,"object")&&"h"in c&&"s"in c&&"l"in c?(b=a.hsl2rgb(c),c.r=b.r,c.g=b.g,c.b=b.b,c.opacity=1,c.hex=b.hex):(y(c,"string")&&(c=a.getRGB(c)),y(c,"object")&&"r"in c&&"g"in c&&"b"in c&&!("error"in c)?(b=a.rgb2hsl(c),c.h=b.h,c.s=b.s,c.l=b.l,b=a.rgb2hsb(c),c.v=b.b):(c={hex:"none"},
c.r=c.g=c.b=c.h=c.s=c.v=c.l=-1,c.error=1));c.toString=ka;return c};a.hsb2rgb=function(c,a,b,m){y(c,"object")&&"h"in c&&"s"in c&&"b"in c&&(b=c.b,a=c.s,c=c.h,m=c.o);var e,h,d;c=360*c%360/60;d=b*a;a=d*(1-Y(c%2-1));b=e=h=b-d;c=~~c;b+=[d,a,0,0,a,d][c];e+=[a,d,d,a,0,0][c];h+=[0,0,a,d,d,a][c];return oa(b,e,h,m)};a.hsl2rgb=function(c,a,b,m){y(c,"object")&&"h"in c&&"s"in c&&"l"in c&&(b=c.l,a=c.s,c=c.h);if(1<c||1<a||1<b)c/=360,a/=100,b/=100;var e,h,d;c=360*c%360/60;d=2*a*(0.5>b?b:1-b);a=d*(1-Y(c%2-1));b=e=
h=b-d/2;c=~~c;b+=[d,a,0,0,a,d][c];e+=[a,d,d,a,0,0][c];h+=[0,0,a,d,d,a][c];return oa(b,e,h,m)};a.rgb2hsb=function(c,a,b){b=D(c,a,b);c=b[0];a=b[1];b=b[2];var m,e;m=P(c,a,b);e=m-Q(c,a,b);c=((0==e?0:m==c?(a-b)/e:m==a?(b-c)/e+2:(c-a)/e+4)+360)%6*60/360;return{h:c,s:0==e?0:e/m,b:m,toString:qa}};a.rgb2hsl=function(c,a,b){b=D(c,a,b);c=b[0];a=b[1];b=b[2];var m,e,h;m=P(c,a,b);e=Q(c,a,b);h=m-e;c=((0==h?0:m==c?(a-b)/h:m==a?(b-c)/h+2:(c-a)/h+4)+360)%6*60/360;m=(m+e)/2;return{h:c,s:0==h?0:0.5>m?h/(2*m):h/(2-2*
m),l:m,toString:ra}};a.parsePathString=function(c){if(!c)return null;var b=a.path(c);if(b.arr)return a.path.clone(b.arr);var m={a:7,c:6,o:2,h:1,l:2,m:2,r:4,q:4,s:4,t:2,v:1,u:3,z:0},e=[];y(c,"array")&&y(c[0],"array")&&(e=a.path.clone(c));e.length||J(c).replace(W,function(c,a,b){var h=[];c=a.toLowerCase();b.replace(Z,function(c,a){a&&h.push(+a)});"m"==c&&2<h.length&&(e.push([a].concat(h.splice(0,2))),c="l",a="m"==a?"l":"L");"o"==c&&1==h.length&&e.push([a,h[0] ]);if("r"==c)e.push([a].concat(h));else for(;h.length>=
m[c]&&(e.push([a].concat(h.splice(0,m[c]))),m[c]););});e.toString=a.path.toString;b.arr=a.path.clone(e);return e};var O=a.parseTransformString=function(c){if(!c)return null;var b=[];y(c,"array")&&y(c[0],"array")&&(b=a.path.clone(c));b.length||J(c).replace(ma,function(c,a,m){var e=[];a.toLowerCase();m.replace(Z,function(c,a){a&&e.push(+a)});b.push([a].concat(e))});b.toString=a.path.toString;return b};a._.svgTransform2string=d;a._.rgTransform=RegExp("^[a-z][\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*-?\\.?\\d",
"i");a._.transform2matrix=f;a._unit2px=b;a._.getSomeDefs=u;a._.getSomeSVG=p;a.select=function(c){return x(G.doc.querySelector(c))};a.selectAll=function(c){c=G.doc.querySelectorAll(c);for(var b=(a.set||Array)(),m=0;m<c.length;m++)b.push(x(c[m]));return b};setInterval(function(){for(var c in E)if(E[h](c)){var a=E[c],b=a.node;("svg"!=a.type&&!b.ownerSVGElement||"svg"==a.type&&(!b.parentNode||"ownerSVGElement"in b.parentNode&&!b.ownerSVGElement))&&delete E[c]}},1E4);(function(c){function m(c){function a(c,
b){var m=v(c.node,b);(m=(m=m&&m.match(d))&&m[2])&&"#"==m.charAt()&&(m=m.substring(1))&&(f[m]=(f[m]||[]).concat(function(a){var m={};m[b]=ca(a);v(c.node,m)}))}function b(c){var a=v(c.node,"xlink:href");a&&"#"==a.charAt()&&(a=a.substring(1))&&(f[a]=(f[a]||[]).concat(function(a){c.attr("xlink:href","#"+a)}))}var e=c.selectAll("*"),h,d=/^\s*url\(("|'|)(.*)\1\)\s*$/;c=[];for(var f={},l=0,E=e.length;l<E;l++){h=e[l];a(h,"fill");a(h,"stroke");a(h,"filter");a(h,"mask");a(h,"clip-path");b(h);var t=v(h.node,
"id");t&&(v(h.node,{id:h.id}),c.push({old:t,id:h.id}))}l=0;for(E=c.length;l<E;l++)if(e=f[c[l].old])for(h=0,t=e.length;h<t;h++)e[h](c[l].id)}function e(c,a,b){return function(m){m=m.slice(c,a);1==m.length&&(m=m[0]);return b?b(m):m}}function d(c){return function(){var a=c?"<"+this.type:"",b=this.node.attributes,m=this.node.childNodes;if(c)for(var e=0,h=b.length;e<h;e++)a+=" "+b[e].name+'="'+b[e].value.replace(/"/g,'\\"')+'"';if(m.length){c&&(a+=">");e=0;for(h=m.length;e<h;e++)3==m[e].nodeType?a+=m[e].nodeValue:
1==m[e].nodeType&&(a+=x(m[e]).toString());c&&(a+="</"+this.type+">")}else c&&(a+="/>");return a}}c.attr=function(c,a){if(!c)return this;if(y(c,"string"))if(1<arguments.length){var b={};b[c]=a;c=b}else return k("snap.util.getattr."+c,this).firstDefined();for(var m in c)c[h](m)&&k("snap.util.attr."+m,this,c[m]);return this};c.getBBox=function(c){if(!a.Matrix||!a.path)return this.node.getBBox();var b=this,m=new a.Matrix;if(b.removed)return a._.box();for(;"use"==b.type;)if(c||(m=m.add(b.transform().localMatrix.translate(b.attr("x")||
0,b.attr("y")||0))),b.original)b=b.original;else var e=b.attr("xlink:href"),b=b.original=b.node.ownerDocument.getElementById(e.substring(e.indexOf("#")+1));var e=b._,h=a.path.get[b.type]||a.path.get.deflt;try{if(c)return e.bboxwt=h?a.path.getBBox(b.realPath=h(b)):a._.box(b.node.getBBox()),a._.box(e.bboxwt);b.realPath=h(b);b.matrix=b.transform().localMatrix;e.bbox=a.path.getBBox(a.path.map(b.realPath,m.add(b.matrix)));return a._.box(e.bbox)}catch(d){return a._.box()}};var f=function(){return this.string};
c.transform=function(c){var b=this._;if(null==c){var m=this;c=new a.Matrix(this.node.getCTM());for(var e=n(this),h=[e],d=new a.Matrix,l=e.toTransformString(),b=J(e)==J(this.matrix)?J(b.transform):l;"svg"!=m.type&&(m=m.parent());)h.push(n(m));for(m=h.length;m--;)d.add(h[m]);return{string:b,globalMatrix:c,totalMatrix:d,localMatrix:e,diffMatrix:c.clone().add(e.invert()),global:c.toTransformString(),total:d.toTransformString(),local:l,toString:f}}c instanceof a.Matrix?this.matrix=c:n(this,c);this.node&&
("linearGradient"==this.type||"radialGradient"==this.type?v(this.node,{gradientTransform:this.matrix}):"pattern"==this.type?v(this.node,{patternTransform:this.matrix}):v(this.node,{transform:this.matrix}));return this};c.parent=function(){return x(this.node.parentNode)};c.append=c.add=function(c){if(c){if("set"==c.type){var a=this;c.forEach(function(c){a.add(c)});return this}c=x(c);this.node.appendChild(c.node);c.paper=this.paper}return this};c.appendTo=function(c){c&&(c=x(c),c.append(this));return this};
c.prepend=function(c){if(c){if("set"==c.type){var a=this,b;c.forEach(function(c){b?b.after(c):a.prepend(c);b=c});return this}c=x(c);var m=c.parent();this.node.insertBefore(c.node,this.node.firstChild);this.add&&this.add();c.paper=this.paper;this.parent()&&this.parent().add();m&&m.add()}return this};c.prependTo=function(c){c=x(c);c.prepend(this);return this};c.before=function(c){if("set"==c.type){var a=this;c.forEach(function(c){var b=c.parent();a.node.parentNode.insertBefore(c.node,a.node);b&&b.add()});
this.parent().add();return this}c=x(c);var b=c.parent();this.node.parentNode.insertBefore(c.node,this.node);this.parent()&&this.parent().add();b&&b.add();c.paper=this.paper;return this};c.after=function(c){c=x(c);var a=c.parent();this.node.nextSibling?this.node.parentNode.insertBefore(c.node,this.node.nextSibling):this.node.parentNode.appendChild(c.node);this.parent()&&this.parent().add();a&&a.add();c.paper=this.paper;return this};c.insertBefore=function(c){c=x(c);var a=this.parent();c.node.parentNode.insertBefore(this.node,
c.node);this.paper=c.paper;a&&a.add();c.parent()&&c.parent().add();return this};c.insertAfter=function(c){c=x(c);var a=this.parent();c.node.parentNode.insertBefore(this.node,c.node.nextSibling);this.paper=c.paper;a&&a.add();c.parent()&&c.parent().add();return this};c.remove=function(){var c=this.parent();this.node.parentNode&&this.node.parentNode.removeChild(this.node);delete this.paper;this.removed=!0;c&&c.add();return this};c.select=function(c){return x(this.node.querySelector(c))};c.selectAll=
function(c){c=this.node.querySelectorAll(c);for(var b=(a.set||Array)(),m=0;m<c.length;m++)b.push(x(c[m]));return b};c.asPX=function(c,a){null==a&&(a=this.attr(c));return+b(this,c,a)};c.use=function(){var c,a=this.node.id;a||(a=this.id,v(this.node,{id:a}));c="linearGradient"==this.type||"radialGradient"==this.type||"pattern"==this.type?r(this.type,this.node.parentNode):r("use",this.node.parentNode);v(c.node,{"xlink:href":"#"+a});c.original=this;return c};var l=/\S+/g;c.addClass=function(c){var a=(c||
"").match(l)||[];c=this.node;var b=c.className.baseVal,m=b.match(l)||[],e,h,d;if(a.length){for(e=0;d=a[e++];)h=m.indexOf(d),~h||m.push(d);a=m.join(" ");b!=a&&(c.className.baseVal=a)}return this};c.removeClass=function(c){var a=(c||"").match(l)||[];c=this.node;var b=c.className.baseVal,m=b.match(l)||[],e,h;if(m.length){for(e=0;h=a[e++];)h=m.indexOf(h),~h&&m.splice(h,1);a=m.join(" ");b!=a&&(c.className.baseVal=a)}return this};c.hasClass=function(c){return!!~(this.node.className.baseVal.match(l)||[]).indexOf(c)};
c.toggleClass=function(c,a){if(null!=a)return a?this.addClass(c):this.removeClass(c);var b=(c||"").match(l)||[],m=this.node,e=m.className.baseVal,h=e.match(l)||[],d,f,E;for(d=0;E=b[d++];)f=h.indexOf(E),~f?h.splice(f,1):h.push(E);b=h.join(" ");e!=b&&(m.className.baseVal=b);return this};c.clone=function(){var c=x(this.node.cloneNode(!0));v(c.node,"id")&&v(c.node,{id:c.id});m(c);c.insertAfter(this);return c};c.toDefs=function(){u(this).appendChild(this.node);return this};c.pattern=c.toPattern=function(c,
a,b,m){var e=r("pattern",u(this));null==c&&(c=this.getBBox());y(c,"object")&&"x"in c&&(a=c.y,b=c.width,m=c.height,c=c.x);v(e.node,{x:c,y:a,width:b,height:m,patternUnits:"userSpaceOnUse",id:e.id,viewBox:[c,a,b,m].join(" ")});e.node.appendChild(this.node);return e};c.marker=function(c,a,b,m,e,h){var d=r("marker",u(this));null==c&&(c=this.getBBox());y(c,"object")&&"x"in c&&(a=c.y,b=c.width,m=c.height,e=c.refX||c.cx,h=c.refY||c.cy,c=c.x);v(d.node,{viewBox:[c,a,b,m].join(" "),markerWidth:b,markerHeight:m,
orient:"auto",refX:e||0,refY:h||0,id:d.id});d.node.appendChild(this.node);return d};var E=function(c,a,b,m){"function"!=typeof b||b.length||(m=b,b=L.linear);this.attr=c;this.dur=a;b&&(this.easing=b);m&&(this.callback=m)};a._.Animation=E;a.animation=function(c,a,b,m){return new E(c,a,b,m)};c.inAnim=function(){var c=[],a;for(a in this.anims)this.anims[h](a)&&function(a){c.push({anim:new E(a._attrs,a.dur,a.easing,a._callback),mina:a,curStatus:a.status(),status:function(c){return a.status(c)},stop:function(){a.stop()}})}(this.anims[a]);
return c};a.animate=function(c,a,b,m,e,h){"function"!=typeof e||e.length||(h=e,e=L.linear);var d=L.time();c=L(c,a,d,d+m,L.time,b,e);h&&k.once("mina.finish."+c.id,h);return c};c.stop=function(){for(var c=this.inAnim(),a=0,b=c.length;a<b;a++)c[a].stop();return this};c.animate=function(c,a,b,m){"function"!=typeof b||b.length||(m=b,b=L.linear);c instanceof E&&(m=c.callback,b=c.easing,a=b.dur,c=c.attr);var d=[],f=[],l={},t,ca,n,T=this,q;for(q in c)if(c[h](q)){T.equal?(n=T.equal(q,J(c[q])),t=n.from,ca=
n.to,n=n.f):(t=+T.attr(q),ca=+c[q]);var la=y(t,"array")?t.length:1;l[q]=e(d.length,d.length+la,n);d=d.concat(t);f=f.concat(ca)}t=L.time();var p=L(d,f,t,t+a,L.time,function(c){var a={},b;for(b in l)l[h](b)&&(a[b]=l[b](c));T.attr(a)},b);T.anims[p.id]=p;p._attrs=c;p._callback=m;k("snap.animcreated."+T.id,p);k.once("mina.finish."+p.id,function(){delete T.anims[p.id];m&&m.call(T)});k.once("mina.stop."+p.id,function(){delete T.anims[p.id]});return T};var T={};c.data=function(c,b){var m=T[this.id]=T[this.id]||
{};if(0==arguments.length)return k("snap.data.get."+this.id,this,m,null),m;if(1==arguments.length){if(a.is(c,"object")){for(var e in c)c[h](e)&&this.data(e,c[e]);return this}k("snap.data.get."+this.id,this,m[c],c);return m[c]}m[c]=b;k("snap.data.set."+this.id,this,b,c);return this};c.removeData=function(c){null==c?T[this.id]={}:T[this.id]&&delete T[this.id][c];return this};c.outerSVG=c.toString=d(1);c.innerSVG=d()})(e.prototype);a.parse=function(c){var a=G.doc.createDocumentFragment(),b=!0,m=G.doc.createElement("div");
c=J(c);c.match(/^\s*<\s*svg(?:\s|>)/)||(c="<svg>"+c+"</svg>",b=!1);m.innerHTML=c;if(c=m.getElementsByTagName("svg")[0])if(b)a=c;else for(;c.firstChild;)a.appendChild(c.firstChild);m.innerHTML=aa;return new l(a)};l.prototype.select=e.prototype.select;l.prototype.selectAll=e.prototype.selectAll;a.fragment=function(){for(var c=Array.prototype.slice.call(arguments,0),b=G.doc.createDocumentFragment(),m=0,e=c.length;m<e;m++){var h=c[m];h.node&&h.node.nodeType&&b.appendChild(h.node);h.nodeType&&b.appendChild(h);
"string"==typeof h&&b.appendChild(a.parse(h).node)}return new l(b)};a._.make=r;a._.wrap=x;s.prototype.el=function(c,a){var b=r(c,this.node);a&&b.attr(a);return b};k.on("snap.util.getattr",function(){var c=k.nt(),c=c.substring(c.lastIndexOf(".")+1),a=c.replace(/[A-Z]/g,function(c){return"-"+c.toLowerCase()});return pa[h](a)?this.node.ownerDocument.defaultView.getComputedStyle(this.node,null).getPropertyValue(a):v(this.node,c)});var pa={"alignment-baseline":0,"baseline-shift":0,clip:0,"clip-path":0,
"clip-rule":0,color:0,"color-interpolation":0,"color-interpolation-filters":0,"color-profile":0,"color-rendering":0,cursor:0,direction:0,display:0,"dominant-baseline":0,"enable-background":0,fill:0,"fill-opacity":0,"fill-rule":0,filter:0,"flood-color":0,"flood-opacity":0,font:0,"font-family":0,"font-size":0,"font-size-adjust":0,"font-stretch":0,"font-style":0,"font-variant":0,"font-weight":0,"glyph-orientation-horizontal":0,"glyph-orientation-vertical":0,"image-rendering":0,kerning:0,"letter-spacing":0,
"lighting-color":0,marker:0,"marker-end":0,"marker-mid":0,"marker-start":0,mask:0,opacity:0,overflow:0,"pointer-events":0,"shape-rendering":0,"stop-color":0,"stop-opacity":0,stroke:0,"stroke-dasharray":0,"stroke-dashoffset":0,"stroke-linecap":0,"stroke-linejoin":0,"stroke-miterlimit":0,"stroke-opacity":0,"stroke-width":0,"text-anchor":0,"text-decoration":0,"text-rendering":0,"unicode-bidi":0,visibility:0,"word-spacing":0,"writing-mode":0};k.on("snap.util.attr",function(c){var a=k.nt(),b={},a=a.substring(a.lastIndexOf(".")+
1);b[a]=c;var m=a.replace(/-(\w)/gi,function(c,a){return a.toUpperCase()}),a=a.replace(/[A-Z]/g,function(c){return"-"+c.toLowerCase()});pa[h](a)?this.node.style[m]=null==c?aa:c:v(this.node,b)});a.ajax=function(c,a,b,m){var e=new XMLHttpRequest,h=V();if(e){if(y(a,"function"))m=b,b=a,a=null;else if(y(a,"object")){var d=[],f;for(f in a)a.hasOwnProperty(f)&&d.push(encodeURIComponent(f)+"="+encodeURIComponent(a[f]));a=d.join("&")}e.open(a?"POST":"GET",c,!0);a&&(e.setRequestHeader("X-Requested-With","XMLHttpRequest"),
e.setRequestHeader("Content-type","application/x-www-form-urlencoded"));b&&(k.once("snap.ajax."+h+".0",b),k.once("snap.ajax."+h+".200",b),k.once("snap.ajax."+h+".304",b));e.onreadystatechange=function(){4==e.readyState&&k("snap.ajax."+h+"."+e.status,m,e)};if(4==e.readyState)return e;e.send(a);return e}};a.load=function(c,b,m){a.ajax(c,function(c){c=a.parse(c.responseText);m?b.call(m,c):b(c)})};a.getElementByPoint=function(c,a){var b,m,e=G.doc.elementFromPoint(c,a);if(G.win.opera&&"svg"==e.tagName){b=
e;m=b.getBoundingClientRect();b=b.ownerDocument;var h=b.body,d=b.documentElement;b=m.top+(g.win.pageYOffset||d.scrollTop||h.scrollTop)-(d.clientTop||h.clientTop||0);m=m.left+(g.win.pageXOffset||d.scrollLeft||h.scrollLeft)-(d.clientLeft||h.clientLeft||0);h=e.createSVGRect();h.x=c-m;h.y=a-b;h.width=h.height=1;b=e.getIntersectionList(h,null);b.length&&(e=b[b.length-1])}return e?x(e):null};a.plugin=function(c){c(a,e,s,G,l)};return G.win.Snap=a}();C.plugin(function(a,k,y,M,A){function w(a,d,f,b,q,e){null==
d&&"[object SVGMatrix]"==z.call(a)?(this.a=a.a,this.b=a.b,this.c=a.c,this.d=a.d,this.e=a.e,this.f=a.f):null!=a?(this.a=+a,this.b=+d,this.c=+f,this.d=+b,this.e=+q,this.f=+e):(this.a=1,this.c=this.b=0,this.d=1,this.f=this.e=0)}var z=Object.prototype.toString,d=String,f=Math;(function(n){function k(a){return a[0]*a[0]+a[1]*a[1]}function p(a){var d=f.sqrt(k(a));a[0]&&(a[0]/=d);a[1]&&(a[1]/=d)}n.add=function(a,d,e,f,n,p){var k=[[],[],[] ],u=[[this.a,this.c,this.e],[this.b,this.d,this.f],[0,0,1] ];d=[[a,
e,n],[d,f,p],[0,0,1] ];a&&a instanceof w&&(d=[[a.a,a.c,a.e],[a.b,a.d,a.f],[0,0,1] ]);for(a=0;3>a;a++)for(e=0;3>e;e++){for(f=n=0;3>f;f++)n+=u[a][f]*d[f][e];k[a][e]=n}this.a=k[0][0];this.b=k[1][0];this.c=k[0][1];this.d=k[1][1];this.e=k[0][2];this.f=k[1][2];return this};n.invert=function(){var a=this.a*this.d-this.b*this.c;return new w(this.d/a,-this.b/a,-this.c/a,this.a/a,(this.c*this.f-this.d*this.e)/a,(this.b*this.e-this.a*this.f)/a)};n.clone=function(){return new w(this.a,this.b,this.c,this.d,this.e,
this.f)};n.translate=function(a,d){return this.add(1,0,0,1,a,d)};n.scale=function(a,d,e,f){null==d&&(d=a);(e||f)&&this.add(1,0,0,1,e,f);this.add(a,0,0,d,0,0);(e||f)&&this.add(1,0,0,1,-e,-f);return this};n.rotate=function(b,d,e){b=a.rad(b);d=d||0;e=e||0;var l=+f.cos(b).toFixed(9);b=+f.sin(b).toFixed(9);this.add(l,b,-b,l,d,e);return this.add(1,0,0,1,-d,-e)};n.x=function(a,d){return a*this.a+d*this.c+this.e};n.y=function(a,d){return a*this.b+d*this.d+this.f};n.get=function(a){return+this[d.fromCharCode(97+
a)].toFixed(4)};n.toString=function(){return"matrix("+[this.get(0),this.get(1),this.get(2),this.get(3),this.get(4),this.get(5)].join()+")"};n.offset=function(){return[this.e.toFixed(4),this.f.toFixed(4)]};n.determinant=function(){return this.a*this.d-this.b*this.c};n.split=function(){var b={};b.dx=this.e;b.dy=this.f;var d=[[this.a,this.c],[this.b,this.d] ];b.scalex=f.sqrt(k(d[0]));p(d[0]);b.shear=d[0][0]*d[1][0]+d[0][1]*d[1][1];d[1]=[d[1][0]-d[0][0]*b.shear,d[1][1]-d[0][1]*b.shear];b.scaley=f.sqrt(k(d[1]));
p(d[1]);b.shear/=b.scaley;0>this.determinant()&&(b.scalex=-b.scalex);var e=-d[0][1],d=d[1][1];0>d?(b.rotate=a.deg(f.acos(d)),0>e&&(b.rotate=360-b.rotate)):b.rotate=a.deg(f.asin(e));b.isSimple=!+b.shear.toFixed(9)&&(b.scalex.toFixed(9)==b.scaley.toFixed(9)||!b.rotate);b.isSuperSimple=!+b.shear.toFixed(9)&&b.scalex.toFixed(9)==b.scaley.toFixed(9)&&!b.rotate;b.noRotation=!+b.shear.toFixed(9)&&!b.rotate;return b};n.toTransformString=function(a){a=a||this.split();if(+a.shear.toFixed(9))return"m"+[this.get(0),
this.get(1),this.get(2),this.get(3),this.get(4),this.get(5)];a.scalex=+a.scalex.toFixed(4);a.scaley=+a.scaley.toFixed(4);a.rotate=+a.rotate.toFixed(4);return(a.dx||a.dy?"t"+[+a.dx.toFixed(4),+a.dy.toFixed(4)]:"")+(1!=a.scalex||1!=a.scaley?"s"+[a.scalex,a.scaley,0,0]:"")+(a.rotate?"r"+[+a.rotate.toFixed(4),0,0]:"")}})(w.prototype);a.Matrix=w;a.matrix=function(a,d,f,b,k,e){return new w(a,d,f,b,k,e)}});C.plugin(function(a,v,y,M,A){function w(h){return function(d){k.stop();d instanceof A&&1==d.node.childNodes.length&&
("radialGradient"==d.node.firstChild.tagName||"linearGradient"==d.node.firstChild.tagName||"pattern"==d.node.firstChild.tagName)&&(d=d.node.firstChild,b(this).appendChild(d),d=u(d));if(d instanceof v)if("radialGradient"==d.type||"linearGradient"==d.type||"pattern"==d.type){d.node.id||e(d.node,{id:d.id});var f=l(d.node.id)}else f=d.attr(h);else f=a.color(d),f.error?(f=a(b(this).ownerSVGElement).gradient(d))?(f.node.id||e(f.node,{id:f.id}),f=l(f.node.id)):f=d:f=r(f);d={};d[h]=f;e(this.node,d);this.node.style[h]=
x}}function z(a){k.stop();a==+a&&(a+="px");this.node.style.fontSize=a}function d(a){var b=[];a=a.childNodes;for(var e=0,f=a.length;e<f;e++){var l=a[e];3==l.nodeType&&b.push(l.nodeValue);"tspan"==l.tagName&&(1==l.childNodes.length&&3==l.firstChild.nodeType?b.push(l.firstChild.nodeValue):b.push(d(l)))}return b}function f(){k.stop();return this.node.style.fontSize}var n=a._.make,u=a._.wrap,p=a.is,b=a._.getSomeDefs,q=/^url\(#?([^)]+)\)$/,e=a._.$,l=a.url,r=String,s=a._.separator,x="";k.on("snap.util.attr.mask",
function(a){if(a instanceof v||a instanceof A){k.stop();a instanceof A&&1==a.node.childNodes.length&&(a=a.node.firstChild,b(this).appendChild(a),a=u(a));if("mask"==a.type)var d=a;else d=n("mask",b(this)),d.node.appendChild(a.node);!d.node.id&&e(d.node,{id:d.id});e(this.node,{mask:l(d.id)})}});(function(a){k.on("snap.util.attr.clip",a);k.on("snap.util.attr.clip-path",a);k.on("snap.util.attr.clipPath",a)})(function(a){if(a instanceof v||a instanceof A){k.stop();if("clipPath"==a.type)var d=a;else d=
n("clipPath",b(this)),d.node.appendChild(a.node),!d.node.id&&e(d.node,{id:d.id});e(this.node,{"clip-path":l(d.id)})}});k.on("snap.util.attr.fill",w("fill"));k.on("snap.util.attr.stroke",w("stroke"));var G=/^([lr])(?:\(([^)]*)\))?(.*)$/i;k.on("snap.util.grad.parse",function(a){a=r(a);var b=a.match(G);if(!b)return null;a=b[1];var e=b[2],b=b[3],e=e.split(/\s*,\s*/).map(function(a){return+a==a?+a:a});1==e.length&&0==e[0]&&(e=[]);b=b.split("-");b=b.map(function(a){a=a.split(":");var b={color:a[0]};a[1]&&
(b.offset=parseFloat(a[1]));return b});return{type:a,params:e,stops:b}});k.on("snap.util.attr.d",function(b){k.stop();p(b,"array")&&p(b[0],"array")&&(b=a.path.toString.call(b));b=r(b);b.match(/[ruo]/i)&&(b=a.path.toAbsolute(b));e(this.node,{d:b})})(-1);k.on("snap.util.attr.#text",function(a){k.stop();a=r(a);for(a=M.doc.createTextNode(a);this.node.firstChild;)this.node.removeChild(this.node.firstChild);this.node.appendChild(a)})(-1);k.on("snap.util.attr.path",function(a){k.stop();this.attr({d:a})})(-1);
k.on("snap.util.attr.class",function(a){k.stop();this.node.className.baseVal=a})(-1);k.on("snap.util.attr.viewBox",function(a){a=p(a,"object")&&"x"in a?[a.x,a.y,a.width,a.height].join(" "):p(a,"array")?a.join(" "):a;e(this.node,{viewBox:a});k.stop()})(-1);k.on("snap.util.attr.transform",function(a){this.transform(a);k.stop()})(-1);k.on("snap.util.attr.r",function(a){"rect"==this.type&&(k.stop(),e(this.node,{rx:a,ry:a}))})(-1);k.on("snap.util.attr.textpath",function(a){k.stop();if("text"==this.type){var d,
f;if(!a&&this.textPath){for(a=this.textPath;a.node.firstChild;)this.node.appendChild(a.node.firstChild);a.remove();delete this.textPath}else if(p(a,"string")?(d=b(this),a=u(d.parentNode).path(a),d.appendChild(a.node),d=a.id,a.attr({id:d})):(a=u(a),a instanceof v&&(d=a.attr("id"),d||(d=a.id,a.attr({id:d})))),d)if(a=this.textPath,f=this.node,a)a.attr({"xlink:href":"#"+d});else{for(a=e("textPath",{"xlink:href":"#"+d});f.firstChild;)a.appendChild(f.firstChild);f.appendChild(a);this.textPath=u(a)}}})(-1);
k.on("snap.util.attr.text",function(a){if("text"==this.type){for(var b=this.node,d=function(a){var b=e("tspan");if(p(a,"array"))for(var f=0;f<a.length;f++)b.appendChild(d(a[f]));else b.appendChild(M.doc.createTextNode(a));b.normalize&&b.normalize();return b};b.firstChild;)b.removeChild(b.firstChild);for(a=d(a);a.firstChild;)b.appendChild(a.firstChild)}k.stop()})(-1);k.on("snap.util.attr.fontSize",z)(-1);k.on("snap.util.attr.font-size",z)(-1);k.on("snap.util.getattr.transform",function(){k.stop();
return this.transform()})(-1);k.on("snap.util.getattr.textpath",function(){k.stop();return this.textPath})(-1);(function(){function b(d){return function(){k.stop();var b=M.doc.defaultView.getComputedStyle(this.node,null).getPropertyValue("marker-"+d);return"none"==b?b:a(M.doc.getElementById(b.match(q)[1]))}}function d(a){return function(b){k.stop();var d="marker"+a.charAt(0).toUpperCase()+a.substring(1);if(""==b||!b)this.node.style[d]="none";else if("marker"==b.type){var f=b.node.id;f||e(b.node,{id:b.id});
this.node.style[d]=l(f)}}}k.on("snap.util.getattr.marker-end",b("end"))(-1);k.on("snap.util.getattr.markerEnd",b("end"))(-1);k.on("snap.util.getattr.marker-start",b("start"))(-1);k.on("snap.util.getattr.markerStart",b("start"))(-1);k.on("snap.util.getattr.marker-mid",b("mid"))(-1);k.on("snap.util.getattr.markerMid",b("mid"))(-1);k.on("snap.util.attr.marker-end",d("end"))(-1);k.on("snap.util.attr.markerEnd",d("end"))(-1);k.on("snap.util.attr.marker-start",d("start"))(-1);k.on("snap.util.attr.markerStart",
d("start"))(-1);k.on("snap.util.attr.marker-mid",d("mid"))(-1);k.on("snap.util.attr.markerMid",d("mid"))(-1)})();k.on("snap.util.getattr.r",function(){if("rect"==this.type&&e(this.node,"rx")==e(this.node,"ry"))return k.stop(),e(this.node,"rx")})(-1);k.on("snap.util.getattr.text",function(){if("text"==this.type||"tspan"==this.type){k.stop();var a=d(this.node);return 1==a.length?a[0]:a}})(-1);k.on("snap.util.getattr.#text",function(){return this.node.textContent})(-1);k.on("snap.util.getattr.viewBox",
function(){k.stop();var b=e(this.node,"viewBox");if(b)return b=b.split(s),a._.box(+b[0],+b[1],+b[2],+b[3])})(-1);k.on("snap.util.getattr.points",function(){var a=e(this.node,"points");k.stop();if(a)return a.split(s)})(-1);k.on("snap.util.getattr.path",function(){var a=e(this.node,"d");k.stop();return a})(-1);k.on("snap.util.getattr.class",function(){return this.node.className.baseVal})(-1);k.on("snap.util.getattr.fontSize",f)(-1);k.on("snap.util.getattr.font-size",f)(-1)});C.plugin(function(a,v,y,
M,A){function w(a){return a}function z(a){return function(b){return+b.toFixed(3)+a}}var d={"+":function(a,b){return a+b},"-":function(a,b){return a-b},"/":function(a,b){return a/b},"*":function(a,b){return a*b}},f=String,n=/[a-z]+$/i,u=/^\s*([+\-\/*])\s*=\s*([\d.eE+\-]+)\s*([^\d\s]+)?\s*$/;k.on("snap.util.attr",function(a){if(a=f(a).match(u)){var b=k.nt(),b=b.substring(b.lastIndexOf(".")+1),q=this.attr(b),e={};k.stop();var l=a[3]||"",r=q.match(n),s=d[a[1] ];r&&r==l?a=s(parseFloat(q),+a[2]):(q=this.asPX(b),
a=s(this.asPX(b),this.asPX(b,a[2]+l)));isNaN(q)||isNaN(a)||(e[b]=a,this.attr(e))}})(-10);k.on("snap.util.equal",function(a,b){var q=f(this.attr(a)||""),e=f(b).match(u);if(e){k.stop();var l=e[3]||"",r=q.match(n),s=d[e[1] ];if(r&&r==l)return{from:parseFloat(q),to:s(parseFloat(q),+e[2]),f:z(r)};q=this.asPX(a);return{from:q,to:s(q,this.asPX(a,e[2]+l)),f:w}}})(-10)});C.plugin(function(a,v,y,M,A){var w=y.prototype,z=a.is;w.rect=function(a,d,k,p,b,q){var e;null==q&&(q=b);z(a,"object")&&"[object Object]"==
a?e=a:null!=a&&(e={x:a,y:d,width:k,height:p},null!=b&&(e.rx=b,e.ry=q));return this.el("rect",e)};w.circle=function(a,d,k){var p;z(a,"object")&&"[object Object]"==a?p=a:null!=a&&(p={cx:a,cy:d,r:k});return this.el("circle",p)};var d=function(){function a(){this.parentNode.removeChild(this)}return function(d,k){var p=M.doc.createElement("img"),b=M.doc.body;p.style.cssText="position:absolute;left:-9999em;top:-9999em";p.onload=function(){k.call(p);p.onload=p.onerror=null;b.removeChild(p)};p.onerror=a;
b.appendChild(p);p.src=d}}();w.image=function(f,n,k,p,b){var q=this.el("image");if(z(f,"object")&&"src"in f)q.attr(f);else if(null!=f){var e={"xlink:href":f,preserveAspectRatio:"none"};null!=n&&null!=k&&(e.x=n,e.y=k);null!=p&&null!=b?(e.width=p,e.height=b):d(f,function(){a._.$(q.node,{width:this.offsetWidth,height:this.offsetHeight})});a._.$(q.node,e)}return q};w.ellipse=function(a,d,k,p){var b;z(a,"object")&&"[object Object]"==a?b=a:null!=a&&(b={cx:a,cy:d,rx:k,ry:p});return this.el("ellipse",b)};
w.path=function(a){var d;z(a,"object")&&!z(a,"array")?d=a:a&&(d={d:a});return this.el("path",d)};w.group=w.g=function(a){var d=this.el("g");1==arguments.length&&a&&!a.type?d.attr(a):arguments.length&&d.add(Array.prototype.slice.call(arguments,0));return d};w.svg=function(a,d,k,p,b,q,e,l){var r={};z(a,"object")&&null==d?r=a:(null!=a&&(r.x=a),null!=d&&(r.y=d),null!=k&&(r.width=k),null!=p&&(r.height=p),null!=b&&null!=q&&null!=e&&null!=l&&(r.viewBox=[b,q,e,l]));return this.el("svg",r)};w.mask=function(a){var d=
this.el("mask");1==arguments.length&&a&&!a.type?d.attr(a):arguments.length&&d.add(Array.prototype.slice.call(arguments,0));return d};w.ptrn=function(a,d,k,p,b,q,e,l){if(z(a,"object"))var r=a;else arguments.length?(r={},null!=a&&(r.x=a),null!=d&&(r.y=d),null!=k&&(r.width=k),null!=p&&(r.height=p),null!=b&&null!=q&&null!=e&&null!=l&&(r.viewBox=[b,q,e,l])):r={patternUnits:"userSpaceOnUse"};return this.el("pattern",r)};w.use=function(a){return null!=a?(make("use",this.node),a instanceof v&&(a.attr("id")||
a.attr({id:ID()}),a=a.attr("id")),this.el("use",{"xlink:href":a})):v.prototype.use.call(this)};w.text=function(a,d,k){var p={};z(a,"object")?p=a:null!=a&&(p={x:a,y:d,text:k||""});return this.el("text",p)};w.line=function(a,d,k,p){var b={};z(a,"object")?b=a:null!=a&&(b={x1:a,x2:k,y1:d,y2:p});return this.el("line",b)};w.polyline=function(a){1<arguments.length&&(a=Array.prototype.slice.call(arguments,0));var d={};z(a,"object")&&!z(a,"array")?d=a:null!=a&&(d={points:a});return this.el("polyline",d)};
w.polygon=function(a){1<arguments.length&&(a=Array.prototype.slice.call(arguments,0));var d={};z(a,"object")&&!z(a,"array")?d=a:null!=a&&(d={points:a});return this.el("polygon",d)};(function(){function d(){return this.selectAll("stop")}function n(b,d){var f=e("stop"),k={offset:+d+"%"};b=a.color(b);k["stop-color"]=b.hex;1>b.opacity&&(k["stop-opacity"]=b.opacity);e(f,k);this.node.appendChild(f);return this}function u(){if("linearGradient"==this.type){var b=e(this.node,"x1")||0,d=e(this.node,"x2")||
1,f=e(this.node,"y1")||0,k=e(this.node,"y2")||0;return a._.box(b,f,math.abs(d-b),math.abs(k-f))}b=this.node.r||0;return a._.box((this.node.cx||0.5)-b,(this.node.cy||0.5)-b,2*b,2*b)}function p(a,d){function f(a,b){for(var d=(b-u)/(a-w),e=w;e<a;e++)h[e].offset=+(+u+d*(e-w)).toFixed(2);w=a;u=b}var n=k("snap.util.grad.parse",null,d).firstDefined(),p;if(!n)return null;n.params.unshift(a);p="l"==n.type.toLowerCase()?b.apply(0,n.params):q.apply(0,n.params);n.type!=n.type.toLowerCase()&&e(p.node,{gradientUnits:"userSpaceOnUse"});
var h=n.stops,n=h.length,u=0,w=0;n--;for(var v=0;v<n;v++)"offset"in h[v]&&f(v,h[v].offset);h[n].offset=h[n].offset||100;f(n,h[n].offset);for(v=0;v<=n;v++){var y=h[v];p.addStop(y.color,y.offset)}return p}function b(b,k,p,q,w){b=a._.make("linearGradient",b);b.stops=d;b.addStop=n;b.getBBox=u;null!=k&&e(b.node,{x1:k,y1:p,x2:q,y2:w});return b}function q(b,k,p,q,w,h){b=a._.make("radialGradient",b);b.stops=d;b.addStop=n;b.getBBox=u;null!=k&&e(b.node,{cx:k,cy:p,r:q});null!=w&&null!=h&&e(b.node,{fx:w,fy:h});
return b}var e=a._.$;w.gradient=function(a){return p(this.defs,a)};w.gradientLinear=function(a,d,e,f){return b(this.defs,a,d,e,f)};w.gradientRadial=function(a,b,d,e,f){return q(this.defs,a,b,d,e,f)};w.toString=function(){var b=this.node.ownerDocument,d=b.createDocumentFragment(),b=b.createElement("div"),e=this.node.cloneNode(!0);d.appendChild(b);b.appendChild(e);a._.$(e,{xmlns:"http://www.w3.org/2000/svg"});b=b.innerHTML;d.removeChild(d.firstChild);return b};w.clear=function(){for(var a=this.node.firstChild,
b;a;)b=a.nextSibling,"defs"!=a.tagName?a.parentNode.removeChild(a):w.clear.call({node:a}),a=b}})()});C.plugin(function(a,k,y,M){function A(a){var b=A.ps=A.ps||{};b[a]?b[a].sleep=100:b[a]={sleep:100};setTimeout(function(){for(var d in b)b[L](d)&&d!=a&&(b[d].sleep--,!b[d].sleep&&delete b[d])});return b[a]}function w(a,b,d,e){null==a&&(a=b=d=e=0);null==b&&(b=a.y,d=a.width,e=a.height,a=a.x);return{x:a,y:b,width:d,w:d,height:e,h:e,x2:a+d,y2:b+e,cx:a+d/2,cy:b+e/2,r1:F.min(d,e)/2,r2:F.max(d,e)/2,r0:F.sqrt(d*
d+e*e)/2,path:s(a,b,d,e),vb:[a,b,d,e].join(" ")}}function z(){return this.join(",").replace(N,"$1")}function d(a){a=C(a);a.toString=z;return a}function f(a,b,d,h,f,k,l,n,p){if(null==p)return e(a,b,d,h,f,k,l,n);if(0>p||e(a,b,d,h,f,k,l,n)<p)p=void 0;else{var q=0.5,O=1-q,s;for(s=e(a,b,d,h,f,k,l,n,O);0.01<Z(s-p);)q/=2,O+=(s<p?1:-1)*q,s=e(a,b,d,h,f,k,l,n,O);p=O}return u(a,b,d,h,f,k,l,n,p)}function n(b,d){function e(a){return+(+a).toFixed(3)}return a._.cacher(function(a,h,l){a instanceof k&&(a=a.attr("d"));
a=I(a);for(var n,p,D,q,O="",s={},c=0,t=0,r=a.length;t<r;t++){D=a[t];if("M"==D[0])n=+D[1],p=+D[2];else{q=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6]);if(c+q>h){if(d&&!s.start){n=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6],h-c);O+=["C"+e(n.start.x),e(n.start.y),e(n.m.x),e(n.m.y),e(n.x),e(n.y)];if(l)return O;s.start=O;O=["M"+e(n.x),e(n.y)+"C"+e(n.n.x),e(n.n.y),e(n.end.x),e(n.end.y),e(D[5]),e(D[6])].join();c+=q;n=+D[5];p=+D[6];continue}if(!b&&!d)return n=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6],h-c)}c+=q;n=+D[5];p=+D[6]}O+=
D.shift()+D}s.end=O;return n=b?c:d?s:u(n,p,D[0],D[1],D[2],D[3],D[4],D[5],1)},null,a._.clone)}function u(a,b,d,e,h,f,k,l,n){var p=1-n,q=ma(p,3),s=ma(p,2),c=n*n,t=c*n,r=q*a+3*s*n*d+3*p*n*n*h+t*k,q=q*b+3*s*n*e+3*p*n*n*f+t*l,s=a+2*n*(d-a)+c*(h-2*d+a),t=b+2*n*(e-b)+c*(f-2*e+b),x=d+2*n*(h-d)+c*(k-2*h+d),c=e+2*n*(f-e)+c*(l-2*f+e);a=p*a+n*d;b=p*b+n*e;h=p*h+n*k;f=p*f+n*l;l=90-180*F.atan2(s-x,t-c)/S;return{x:r,y:q,m:{x:s,y:t},n:{x:x,y:c},start:{x:a,y:b},end:{x:h,y:f},alpha:l}}function p(b,d,e,h,f,n,k,l){a.is(b,
"array")||(b=[b,d,e,h,f,n,k,l]);b=U.apply(null,b);return w(b.min.x,b.min.y,b.max.x-b.min.x,b.max.y-b.min.y)}function b(a,b,d){return b>=a.x&&b<=a.x+a.width&&d>=a.y&&d<=a.y+a.height}function q(a,d){a=w(a);d=w(d);return b(d,a.x,a.y)||b(d,a.x2,a.y)||b(d,a.x,a.y2)||b(d,a.x2,a.y2)||b(a,d.x,d.y)||b(a,d.x2,d.y)||b(a,d.x,d.y2)||b(a,d.x2,d.y2)||(a.x<d.x2&&a.x>d.x||d.x<a.x2&&d.x>a.x)&&(a.y<d.y2&&a.y>d.y||d.y<a.y2&&d.y>a.y)}function e(a,b,d,e,h,f,n,k,l){null==l&&(l=1);l=(1<l?1:0>l?0:l)/2;for(var p=[-0.1252,
0.1252,-0.3678,0.3678,-0.5873,0.5873,-0.7699,0.7699,-0.9041,0.9041,-0.9816,0.9816],q=[0.2491,0.2491,0.2335,0.2335,0.2032,0.2032,0.1601,0.1601,0.1069,0.1069,0.0472,0.0472],s=0,c=0;12>c;c++)var t=l*p[c]+l,r=t*(t*(-3*a+9*d-9*h+3*n)+6*a-12*d+6*h)-3*a+3*d,t=t*(t*(-3*b+9*e-9*f+3*k)+6*b-12*e+6*f)-3*b+3*e,s=s+q[c]*F.sqrt(r*r+t*t);return l*s}function l(a,b,d){a=I(a);b=I(b);for(var h,f,l,n,k,s,r,O,x,c,t=d?0:[],w=0,v=a.length;w<v;w++)if(x=a[w],"M"==x[0])h=k=x[1],f=s=x[2];else{"C"==x[0]?(x=[h,f].concat(x.slice(1)),
h=x[6],f=x[7]):(x=[h,f,h,f,k,s,k,s],h=k,f=s);for(var G=0,y=b.length;G<y;G++)if(c=b[G],"M"==c[0])l=r=c[1],n=O=c[2];else{"C"==c[0]?(c=[l,n].concat(c.slice(1)),l=c[6],n=c[7]):(c=[l,n,l,n,r,O,r,O],l=r,n=O);var z;var K=x,B=c;z=d;var H=p(K),J=p(B);if(q(H,J)){for(var H=e.apply(0,K),J=e.apply(0,B),H=~~(H/8),J=~~(J/8),U=[],A=[],F={},M=z?0:[],P=0;P<H+1;P++){var C=u.apply(0,K.concat(P/H));U.push({x:C.x,y:C.y,t:P/H})}for(P=0;P<J+1;P++)C=u.apply(0,B.concat(P/J)),A.push({x:C.x,y:C.y,t:P/J});for(P=0;P<H;P++)for(K=
0;K<J;K++){var Q=U[P],L=U[P+1],B=A[K],C=A[K+1],N=0.001>Z(L.x-Q.x)?"y":"x",S=0.001>Z(C.x-B.x)?"y":"x",R;R=Q.x;var Y=Q.y,V=L.x,ea=L.y,fa=B.x,ga=B.y,ha=C.x,ia=C.y;if(W(R,V)<X(fa,ha)||X(R,V)>W(fa,ha)||W(Y,ea)<X(ga,ia)||X(Y,ea)>W(ga,ia))R=void 0;else{var $=(R*ea-Y*V)*(fa-ha)-(R-V)*(fa*ia-ga*ha),aa=(R*ea-Y*V)*(ga-ia)-(Y-ea)*(fa*ia-ga*ha),ja=(R-V)*(ga-ia)-(Y-ea)*(fa-ha);if(ja){var $=$/ja,aa=aa/ja,ja=+$.toFixed(2),ba=+aa.toFixed(2);R=ja<+X(R,V).toFixed(2)||ja>+W(R,V).toFixed(2)||ja<+X(fa,ha).toFixed(2)||
ja>+W(fa,ha).toFixed(2)||ba<+X(Y,ea).toFixed(2)||ba>+W(Y,ea).toFixed(2)||ba<+X(ga,ia).toFixed(2)||ba>+W(ga,ia).toFixed(2)?void 0:{x:$,y:aa}}else R=void 0}R&&F[R.x.toFixed(4)]!=R.y.toFixed(4)&&(F[R.x.toFixed(4)]=R.y.toFixed(4),Q=Q.t+Z((R[N]-Q[N])/(L[N]-Q[N]))*(L.t-Q.t),B=B.t+Z((R[S]-B[S])/(C[S]-B[S]))*(C.t-B.t),0<=Q&&1>=Q&&0<=B&&1>=B&&(z?M++:M.push({x:R.x,y:R.y,t1:Q,t2:B})))}z=M}else z=z?0:[];if(d)t+=z;else{H=0;for(J=z.length;H<J;H++)z[H].segment1=w,z[H].segment2=G,z[H].bez1=x,z[H].bez2=c;t=t.concat(z)}}}return t}
function r(a){var b=A(a);if(b.bbox)return C(b.bbox);if(!a)return w();a=I(a);for(var d=0,e=0,h=[],f=[],l,n=0,k=a.length;n<k;n++)l=a[n],"M"==l[0]?(d=l[1],e=l[2],h.push(d),f.push(e)):(d=U(d,e,l[1],l[2],l[3],l[4],l[5],l[6]),h=h.concat(d.min.x,d.max.x),f=f.concat(d.min.y,d.max.y),d=l[5],e=l[6]);a=X.apply(0,h);l=X.apply(0,f);h=W.apply(0,h);f=W.apply(0,f);f=w(a,l,h-a,f-l);b.bbox=C(f);return f}function s(a,b,d,e,h){if(h)return[["M",+a+ +h,b],["l",d-2*h,0],["a",h,h,0,0,1,h,h],["l",0,e-2*h],["a",h,h,0,0,1,
-h,h],["l",2*h-d,0],["a",h,h,0,0,1,-h,-h],["l",0,2*h-e],["a",h,h,0,0,1,h,-h],["z"] ];a=[["M",a,b],["l",d,0],["l",0,e],["l",-d,0],["z"] ];a.toString=z;return a}function x(a,b,d,e,h){null==h&&null==e&&(e=d);a=+a;b=+b;d=+d;e=+e;if(null!=h){var f=Math.PI/180,l=a+d*Math.cos(-e*f);a+=d*Math.cos(-h*f);var n=b+d*Math.sin(-e*f);b+=d*Math.sin(-h*f);d=[["M",l,n],["A",d,d,0,+(180<h-e),0,a,b] ]}else d=[["M",a,b],["m",0,-e],["a",d,e,0,1,1,0,2*e],["a",d,e,0,1,1,0,-2*e],["z"] ];d.toString=z;return d}function G(b){var e=
A(b);if(e.abs)return d(e.abs);Q(b,"array")&&Q(b&&b[0],"array")||(b=a.parsePathString(b));if(!b||!b.length)return[["M",0,0] ];var h=[],f=0,l=0,n=0,k=0,p=0;"M"==b[0][0]&&(f=+b[0][1],l=+b[0][2],n=f,k=l,p++,h[0]=["M",f,l]);for(var q=3==b.length&&"M"==b[0][0]&&"R"==b[1][0].toUpperCase()&&"Z"==b[2][0].toUpperCase(),s,r,w=p,c=b.length;w<c;w++){h.push(s=[]);r=b[w];p=r[0];if(p!=p.toUpperCase())switch(s[0]=p.toUpperCase(),s[0]){case "A":s[1]=r[1];s[2]=r[2];s[3]=r[3];s[4]=r[4];s[5]=r[5];s[6]=+r[6]+f;s[7]=+r[7]+
l;break;case "V":s[1]=+r[1]+l;break;case "H":s[1]=+r[1]+f;break;case "R":for(var t=[f,l].concat(r.slice(1)),u=2,v=t.length;u<v;u++)t[u]=+t[u]+f,t[++u]=+t[u]+l;h.pop();h=h.concat(P(t,q));break;case "O":h.pop();t=x(f,l,r[1],r[2]);t.push(t[0]);h=h.concat(t);break;case "U":h.pop();h=h.concat(x(f,l,r[1],r[2],r[3]));s=["U"].concat(h[h.length-1].slice(-2));break;case "M":n=+r[1]+f,k=+r[2]+l;default:for(u=1,v=r.length;u<v;u++)s[u]=+r[u]+(u%2?f:l)}else if("R"==p)t=[f,l].concat(r.slice(1)),h.pop(),h=h.concat(P(t,
q)),s=["R"].concat(r.slice(-2));else if("O"==p)h.pop(),t=x(f,l,r[1],r[2]),t.push(t[0]),h=h.concat(t);else if("U"==p)h.pop(),h=h.concat(x(f,l,r[1],r[2],r[3])),s=["U"].concat(h[h.length-1].slice(-2));else for(t=0,u=r.length;t<u;t++)s[t]=r[t];p=p.toUpperCase();if("O"!=p)switch(s[0]){case "Z":f=+n;l=+k;break;case "H":f=s[1];break;case "V":l=s[1];break;case "M":n=s[s.length-2],k=s[s.length-1];default:f=s[s.length-2],l=s[s.length-1]}}h.toString=z;e.abs=d(h);return h}function h(a,b,d,e){return[a,b,d,e,d,
e]}function J(a,b,d,e,h,f){var l=1/3,n=2/3;return[l*a+n*d,l*b+n*e,l*h+n*d,l*f+n*e,h,f]}function K(b,d,e,h,f,l,n,k,p,s){var r=120*S/180,q=S/180*(+f||0),c=[],t,x=a._.cacher(function(a,b,c){var d=a*F.cos(c)-b*F.sin(c);a=a*F.sin(c)+b*F.cos(c);return{x:d,y:a}});if(s)v=s[0],t=s[1],l=s[2],u=s[3];else{t=x(b,d,-q);b=t.x;d=t.y;t=x(k,p,-q);k=t.x;p=t.y;F.cos(S/180*f);F.sin(S/180*f);t=(b-k)/2;v=(d-p)/2;u=t*t/(e*e)+v*v/(h*h);1<u&&(u=F.sqrt(u),e*=u,h*=u);var u=e*e,w=h*h,u=(l==n?-1:1)*F.sqrt(Z((u*w-u*v*v-w*t*t)/
(u*v*v+w*t*t)));l=u*e*v/h+(b+k)/2;var u=u*-h*t/e+(d+p)/2,v=F.asin(((d-u)/h).toFixed(9));t=F.asin(((p-u)/h).toFixed(9));v=b<l?S-v:v;t=k<l?S-t:t;0>v&&(v=2*S+v);0>t&&(t=2*S+t);n&&v>t&&(v-=2*S);!n&&t>v&&(t-=2*S)}if(Z(t-v)>r){var c=t,w=k,G=p;t=v+r*(n&&t>v?1:-1);k=l+e*F.cos(t);p=u+h*F.sin(t);c=K(k,p,e,h,f,0,n,w,G,[t,c,l,u])}l=t-v;f=F.cos(v);r=F.sin(v);n=F.cos(t);t=F.sin(t);l=F.tan(l/4);e=4/3*e*l;l*=4/3*h;h=[b,d];b=[b+e*r,d-l*f];d=[k+e*t,p-l*n];k=[k,p];b[0]=2*h[0]-b[0];b[1]=2*h[1]-b[1];if(s)return[b,d,k].concat(c);
c=[b,d,k].concat(c).join().split(",");s=[];k=0;for(p=c.length;k<p;k++)s[k]=k%2?x(c[k-1],c[k],q).y:x(c[k],c[k+1],q).x;return s}function U(a,b,d,e,h,f,l,k){for(var n=[],p=[[],[] ],s,r,c,t,q=0;2>q;++q)0==q?(r=6*a-12*d+6*h,s=-3*a+9*d-9*h+3*l,c=3*d-3*a):(r=6*b-12*e+6*f,s=-3*b+9*e-9*f+3*k,c=3*e-3*b),1E-12>Z(s)?1E-12>Z(r)||(s=-c/r,0<s&&1>s&&n.push(s)):(t=r*r-4*c*s,c=F.sqrt(t),0>t||(t=(-r+c)/(2*s),0<t&&1>t&&n.push(t),s=(-r-c)/(2*s),0<s&&1>s&&n.push(s)));for(r=q=n.length;q--;)s=n[q],c=1-s,p[0][q]=c*c*c*a+3*
c*c*s*d+3*c*s*s*h+s*s*s*l,p[1][q]=c*c*c*b+3*c*c*s*e+3*c*s*s*f+s*s*s*k;p[0][r]=a;p[1][r]=b;p[0][r+1]=l;p[1][r+1]=k;p[0].length=p[1].length=r+2;return{min:{x:X.apply(0,p[0]),y:X.apply(0,p[1])},max:{x:W.apply(0,p[0]),y:W.apply(0,p[1])}}}function I(a,b){var e=!b&&A(a);if(!b&&e.curve)return d(e.curve);var f=G(a),l=b&&G(b),n={x:0,y:0,bx:0,by:0,X:0,Y:0,qx:null,qy:null},k={x:0,y:0,bx:0,by:0,X:0,Y:0,qx:null,qy:null},p=function(a,b,c){if(!a)return["C",b.x,b.y,b.x,b.y,b.x,b.y];a[0]in{T:1,Q:1}||(b.qx=b.qy=null);
switch(a[0]){case "M":b.X=a[1];b.Y=a[2];break;case "A":a=["C"].concat(K.apply(0,[b.x,b.y].concat(a.slice(1))));break;case "S":"C"==c||"S"==c?(c=2*b.x-b.bx,b=2*b.y-b.by):(c=b.x,b=b.y);a=["C",c,b].concat(a.slice(1));break;case "T":"Q"==c||"T"==c?(b.qx=2*b.x-b.qx,b.qy=2*b.y-b.qy):(b.qx=b.x,b.qy=b.y);a=["C"].concat(J(b.x,b.y,b.qx,b.qy,a[1],a[2]));break;case "Q":b.qx=a[1];b.qy=a[2];a=["C"].concat(J(b.x,b.y,a[1],a[2],a[3],a[4]));break;case "L":a=["C"].concat(h(b.x,b.y,a[1],a[2]));break;case "H":a=["C"].concat(h(b.x,
b.y,a[1],b.y));break;case "V":a=["C"].concat(h(b.x,b.y,b.x,a[1]));break;case "Z":a=["C"].concat(h(b.x,b.y,b.X,b.Y))}return a},s=function(a,b){if(7<a[b].length){a[b].shift();for(var c=a[b];c.length;)q[b]="A",l&&(u[b]="A"),a.splice(b++,0,["C"].concat(c.splice(0,6)));a.splice(b,1);v=W(f.length,l&&l.length||0)}},r=function(a,b,c,d,e){a&&b&&"M"==a[e][0]&&"M"!=b[e][0]&&(b.splice(e,0,["M",d.x,d.y]),c.bx=0,c.by=0,c.x=a[e][1],c.y=a[e][2],v=W(f.length,l&&l.length||0))},q=[],u=[],c="",t="",x=0,v=W(f.length,
l&&l.length||0);for(;x<v;x++){f[x]&&(c=f[x][0]);"C"!=c&&(q[x]=c,x&&(t=q[x-1]));f[x]=p(f[x],n,t);"A"!=q[x]&&"C"==c&&(q[x]="C");s(f,x);l&&(l[x]&&(c=l[x][0]),"C"!=c&&(u[x]=c,x&&(t=u[x-1])),l[x]=p(l[x],k,t),"A"!=u[x]&&"C"==c&&(u[x]="C"),s(l,x));r(f,l,n,k,x);r(l,f,k,n,x);var w=f[x],z=l&&l[x],y=w.length,U=l&&z.length;n.x=w[y-2];n.y=w[y-1];n.bx=$(w[y-4])||n.x;n.by=$(w[y-3])||n.y;k.bx=l&&($(z[U-4])||k.x);k.by=l&&($(z[U-3])||k.y);k.x=l&&z[U-2];k.y=l&&z[U-1]}l||(e.curve=d(f));return l?[f,l]:f}function P(a,
b){for(var d=[],e=0,h=a.length;h-2*!b>e;e+=2){var f=[{x:+a[e-2],y:+a[e-1]},{x:+a[e],y:+a[e+1]},{x:+a[e+2],y:+a[e+3]},{x:+a[e+4],y:+a[e+5]}];b?e?h-4==e?f[3]={x:+a[0],y:+a[1]}:h-2==e&&(f[2]={x:+a[0],y:+a[1]},f[3]={x:+a[2],y:+a[3]}):f[0]={x:+a[h-2],y:+a[h-1]}:h-4==e?f[3]=f[2]:e||(f[0]={x:+a[e],y:+a[e+1]});d.push(["C",(-f[0].x+6*f[1].x+f[2].x)/6,(-f[0].y+6*f[1].y+f[2].y)/6,(f[1].x+6*f[2].x-f[3].x)/6,(f[1].y+6*f[2].y-f[3].y)/6,f[2].x,f[2].y])}return d}y=k.prototype;var Q=a.is,C=a._.clone,L="hasOwnProperty",
N=/,?([a-z]),?/gi,$=parseFloat,F=Math,S=F.PI,X=F.min,W=F.max,ma=F.pow,Z=F.abs;M=n(1);var na=n(),ba=n(0,1),V=a._unit2px;a.path=A;a.path.getTotalLength=M;a.path.getPointAtLength=na;a.path.getSubpath=function(a,b,d){if(1E-6>this.getTotalLength(a)-d)return ba(a,b).end;a=ba(a,d,1);return b?ba(a,b).end:a};y.getTotalLength=function(){if(this.node.getTotalLength)return this.node.getTotalLength()};y.getPointAtLength=function(a){return na(this.attr("d"),a)};y.getSubpath=function(b,d){return a.path.getSubpath(this.attr("d"),
b,d)};a._.box=w;a.path.findDotsAtSegment=u;a.path.bezierBBox=p;a.path.isPointInsideBBox=b;a.path.isBBoxIntersect=q;a.path.intersection=function(a,b){return l(a,b)};a.path.intersectionNumber=function(a,b){return l(a,b,1)};a.path.isPointInside=function(a,d,e){var h=r(a);return b(h,d,e)&&1==l(a,[["M",d,e],["H",h.x2+10] ],1)%2};a.path.getBBox=r;a.path.get={path:function(a){return a.attr("path")},circle:function(a){a=V(a);return x(a.cx,a.cy,a.r)},ellipse:function(a){a=V(a);return x(a.cx||0,a.cy||0,a.rx,
a.ry)},rect:function(a){a=V(a);return s(a.x||0,a.y||0,a.width,a.height,a.rx,a.ry)},image:function(a){a=V(a);return s(a.x||0,a.y||0,a.width,a.height)},line:function(a){return"M"+[a.attr("x1")||0,a.attr("y1")||0,a.attr("x2"),a.attr("y2")]},polyline:function(a){return"M"+a.attr("points")},polygon:function(a){return"M"+a.attr("points")+"z"},deflt:function(a){a=a.node.getBBox();return s(a.x,a.y,a.width,a.height)}};a.path.toRelative=function(b){var e=A(b),h=String.prototype.toLowerCase;if(e.rel)return d(e.rel);
a.is(b,"array")&&a.is(b&&b[0],"array")||(b=a.parsePathString(b));var f=[],l=0,n=0,k=0,p=0,s=0;"M"==b[0][0]&&(l=b[0][1],n=b[0][2],k=l,p=n,s++,f.push(["M",l,n]));for(var r=b.length;s<r;s++){var q=f[s]=[],x=b[s];if(x[0]!=h.call(x[0]))switch(q[0]=h.call(x[0]),q[0]){case "a":q[1]=x[1];q[2]=x[2];q[3]=x[3];q[4]=x[4];q[5]=x[5];q[6]=+(x[6]-l).toFixed(3);q[7]=+(x[7]-n).toFixed(3);break;case "v":q[1]=+(x[1]-n).toFixed(3);break;case "m":k=x[1],p=x[2];default:for(var c=1,t=x.length;c<t;c++)q[c]=+(x[c]-(c%2?l:
n)).toFixed(3)}else for(f[s]=[],"m"==x[0]&&(k=x[1]+l,p=x[2]+n),q=0,c=x.length;q<c;q++)f[s][q]=x[q];x=f[s].length;switch(f[s][0]){case "z":l=k;n=p;break;case "h":l+=+f[s][x-1];break;case "v":n+=+f[s][x-1];break;default:l+=+f[s][x-2],n+=+f[s][x-1]}}f.toString=z;e.rel=d(f);return f};a.path.toAbsolute=G;a.path.toCubic=I;a.path.map=function(a,b){if(!b)return a;var d,e,h,f,l,n,k;a=I(a);h=0;for(l=a.length;h<l;h++)for(k=a[h],f=1,n=k.length;f<n;f+=2)d=b.x(k[f],k[f+1]),e=b.y(k[f],k[f+1]),k[f]=d,k[f+1]=e;return a};
a.path.toString=z;a.path.clone=d});C.plugin(function(a,v,y,C){var A=Math.max,w=Math.min,z=function(a){this.items=[];this.bindings={};this.length=0;this.type="set";if(a)for(var f=0,n=a.length;f<n;f++)a[f]&&(this[this.items.length]=this.items[this.items.length]=a[f],this.length++)};v=z.prototype;v.push=function(){for(var a,f,n=0,k=arguments.length;n<k;n++)if(a=arguments[n])f=this.items.length,this[f]=this.items[f]=a,this.length++;return this};v.pop=function(){this.length&&delete this[this.length--];
return this.items.pop()};v.forEach=function(a,f){for(var n=0,k=this.items.length;n<k&&!1!==a.call(f,this.items[n],n);n++);return this};v.animate=function(d,f,n,u){"function"!=typeof n||n.length||(u=n,n=L.linear);d instanceof a._.Animation&&(u=d.callback,n=d.easing,f=n.dur,d=d.attr);var p=arguments;if(a.is(d,"array")&&a.is(p[p.length-1],"array"))var b=!0;var q,e=function(){q?this.b=q:q=this.b},l=0,r=u&&function(){l++==this.length&&u.call(this)};return this.forEach(function(a,l){k.once("snap.animcreated."+
a.id,e);b?p[l]&&a.animate.apply(a,p[l]):a.animate(d,f,n,r)})};v.remove=function(){for(;this.length;)this.pop().remove();return this};v.bind=function(a,f,k){var u={};if("function"==typeof f)this.bindings[a]=f;else{var p=k||a;this.bindings[a]=function(a){u[p]=a;f.attr(u)}}return this};v.attr=function(a){var f={},k;for(k in a)if(this.bindings[k])this.bindings[k](a[k]);else f[k]=a[k];a=0;for(k=this.items.length;a<k;a++)this.items[a].attr(f);return this};v.clear=function(){for(;this.length;)this.pop()};
v.splice=function(a,f,k){a=0>a?A(this.length+a,0):a;f=A(0,w(this.length-a,f));var u=[],p=[],b=[],q;for(q=2;q<arguments.length;q++)b.push(arguments[q]);for(q=0;q<f;q++)p.push(this[a+q]);for(;q<this.length-a;q++)u.push(this[a+q]);var e=b.length;for(q=0;q<e+u.length;q++)this.items[a+q]=this[a+q]=q<e?b[q]:u[q-e];for(q=this.items.length=this.length-=f-e;this[q];)delete this[q++];return new z(p)};v.exclude=function(a){for(var f=0,k=this.length;f<k;f++)if(this[f]==a)return this.splice(f,1),!0;return!1};
v.insertAfter=function(a){for(var f=this.items.length;f--;)this.items[f].insertAfter(a);return this};v.getBBox=function(){for(var a=[],f=[],k=[],u=[],p=this.items.length;p--;)if(!this.items[p].removed){var b=this.items[p].getBBox();a.push(b.x);f.push(b.y);k.push(b.x+b.width);u.push(b.y+b.height)}a=w.apply(0,a);f=w.apply(0,f);k=A.apply(0,k);u=A.apply(0,u);return{x:a,y:f,x2:k,y2:u,width:k-a,height:u-f,cx:a+(k-a)/2,cy:f+(u-f)/2}};v.clone=function(a){a=new z;for(var f=0,k=this.items.length;f<k;f++)a.push(this.items[f].clone());
return a};v.toString=function(){return"Snap\u2018s set"};v.type="set";a.set=function(){var a=new z;arguments.length&&a.push.apply(a,Array.prototype.slice.call(arguments,0));return a}});C.plugin(function(a,v,y,C){function A(a){var b=a[0];switch(b.toLowerCase()){case "t":return[b,0,0];case "m":return[b,1,0,0,1,0,0];case "r":return 4==a.length?[b,0,a[2],a[3] ]:[b,0];case "s":return 5==a.length?[b,1,1,a[3],a[4] ]:3==a.length?[b,1,1]:[b,1]}}function w(b,d,f){d=q(d).replace(/\.{3}|\u2026/g,b);b=a.parseTransformString(b)||
[];d=a.parseTransformString(d)||[];for(var k=Math.max(b.length,d.length),p=[],v=[],h=0,w,z,y,I;h<k;h++){y=b[h]||A(d[h]);I=d[h]||A(y);if(y[0]!=I[0]||"r"==y[0].toLowerCase()&&(y[2]!=I[2]||y[3]!=I[3])||"s"==y[0].toLowerCase()&&(y[3]!=I[3]||y[4]!=I[4])){b=a._.transform2matrix(b,f());d=a._.transform2matrix(d,f());p=[["m",b.a,b.b,b.c,b.d,b.e,b.f] ];v=[["m",d.a,d.b,d.c,d.d,d.e,d.f] ];break}p[h]=[];v[h]=[];w=0;for(z=Math.max(y.length,I.length);w<z;w++)w in y&&(p[h][w]=y[w]),w in I&&(v[h][w]=I[w])}return{from:u(p),
to:u(v),f:n(p)}}function z(a){return a}function d(a){return function(b){return+b.toFixed(3)+a}}function f(b){return a.rgb(b[0],b[1],b[2])}function n(a){var b=0,d,f,k,n,h,p,q=[];d=0;for(f=a.length;d<f;d++){h="[";p=['"'+a[d][0]+'"'];k=1;for(n=a[d].length;k<n;k++)p[k]="val["+b++ +"]";h+=p+"]";q[d]=h}return Function("val","return Snap.path.toString.call(["+q+"])")}function u(a){for(var b=[],d=0,f=a.length;d<f;d++)for(var k=1,n=a[d].length;k<n;k++)b.push(a[d][k]);return b}var p={},b=/[a-z]+$/i,q=String;
p.stroke=p.fill="colour";v.prototype.equal=function(a,b){return k("snap.util.equal",this,a,b).firstDefined()};k.on("snap.util.equal",function(e,k){var r,s;r=q(this.attr(e)||"");var x=this;if(r==+r&&k==+k)return{from:+r,to:+k,f:z};if("colour"==p[e])return r=a.color(r),s=a.color(k),{from:[r.r,r.g,r.b,r.opacity],to:[s.r,s.g,s.b,s.opacity],f:f};if("transform"==e||"gradientTransform"==e||"patternTransform"==e)return k instanceof a.Matrix&&(k=k.toTransformString()),a._.rgTransform.test(k)||(k=a._.svgTransform2string(k)),
w(r,k,function(){return x.getBBox(1)});if("d"==e||"path"==e)return r=a.path.toCubic(r,k),{from:u(r[0]),to:u(r[1]),f:n(r[0])};if("points"==e)return r=q(r).split(a._.separator),s=q(k).split(a._.separator),{from:r,to:s,f:function(a){return a}};aUnit=r.match(b);s=q(k).match(b);return aUnit&&aUnit==s?{from:parseFloat(r),to:parseFloat(k),f:d(aUnit)}:{from:this.asPX(e),to:this.asPX(e,k),f:z}})});C.plugin(function(a,v,y,C){var A=v.prototype,w="createTouch"in C.doc;v="click dblclick mousedown mousemove mouseout mouseover mouseup touchstart touchmove touchend touchcancel".split(" ");
var z={mousedown:"touchstart",mousemove:"touchmove",mouseup:"touchend"},d=function(a,b){var d="y"==a?"scrollTop":"scrollLeft",e=b&&b.node?b.node.ownerDocument:C.doc;return e[d in e.documentElement?"documentElement":"body"][d]},f=function(){this.returnValue=!1},n=function(){return this.originalEvent.preventDefault()},u=function(){this.cancelBubble=!0},p=function(){return this.originalEvent.stopPropagation()},b=function(){if(C.doc.addEventListener)return function(a,b,e,f){var k=w&&z[b]?z[b]:b,l=function(k){var l=
d("y",f),q=d("x",f);if(w&&z.hasOwnProperty(b))for(var r=0,u=k.targetTouches&&k.targetTouches.length;r<u;r++)if(k.targetTouches[r].target==a||a.contains(k.targetTouches[r].target)){u=k;k=k.targetTouches[r];k.originalEvent=u;k.preventDefault=n;k.stopPropagation=p;break}return e.call(f,k,k.clientX+q,k.clientY+l)};b!==k&&a.addEventListener(b,l,!1);a.addEventListener(k,l,!1);return function(){b!==k&&a.removeEventListener(b,l,!1);a.removeEventListener(k,l,!1);return!0}};if(C.doc.attachEvent)return function(a,
b,e,h){var k=function(a){a=a||h.node.ownerDocument.window.event;var b=d("y",h),k=d("x",h),k=a.clientX+k,b=a.clientY+b;a.preventDefault=a.preventDefault||f;a.stopPropagation=a.stopPropagation||u;return e.call(h,a,k,b)};a.attachEvent("on"+b,k);return function(){a.detachEvent("on"+b,k);return!0}}}(),q=[],e=function(a){for(var b=a.clientX,e=a.clientY,f=d("y"),l=d("x"),n,p=q.length;p--;){n=q[p];if(w)for(var r=a.touches&&a.touches.length,u;r--;){if(u=a.touches[r],u.identifier==n.el._drag.id||n.el.node.contains(u.target)){b=
u.clientX;e=u.clientY;(a.originalEvent?a.originalEvent:a).preventDefault();break}}else a.preventDefault();b+=l;e+=f;k("snap.drag.move."+n.el.id,n.move_scope||n.el,b-n.el._drag.x,e-n.el._drag.y,b,e,a)}},l=function(b){a.unmousemove(e).unmouseup(l);for(var d=q.length,f;d--;)f=q[d],f.el._drag={},k("snap.drag.end."+f.el.id,f.end_scope||f.start_scope||f.move_scope||f.el,b);q=[]};for(y=v.length;y--;)(function(d){a[d]=A[d]=function(e,f){a.is(e,"function")&&(this.events=this.events||[],this.events.push({name:d,
f:e,unbind:b(this.node||document,d,e,f||this)}));return this};a["un"+d]=A["un"+d]=function(a){for(var b=this.events||[],e=b.length;e--;)if(b[e].name==d&&(b[e].f==a||!a)){b[e].unbind();b.splice(e,1);!b.length&&delete this.events;break}return this}})(v[y]);A.hover=function(a,b,d,e){return this.mouseover(a,d).mouseout(b,e||d)};A.unhover=function(a,b){return this.unmouseover(a).unmouseout(b)};var r=[];A.drag=function(b,d,f,h,n,p){function u(r,v,w){(r.originalEvent||r).preventDefault();this._drag.x=v;
this._drag.y=w;this._drag.id=r.identifier;!q.length&&a.mousemove(e).mouseup(l);q.push({el:this,move_scope:h,start_scope:n,end_scope:p});d&&k.on("snap.drag.start."+this.id,d);b&&k.on("snap.drag.move."+this.id,b);f&&k.on("snap.drag.end."+this.id,f);k("snap.drag.start."+this.id,n||h||this,v,w,r)}if(!arguments.length){var v;return this.drag(function(a,b){this.attr({transform:v+(v?"T":"t")+[a,b]})},function(){v=this.transform().local})}this._drag={};r.push({el:this,start:u});this.mousedown(u);return this};
A.undrag=function(){for(var b=r.length;b--;)r[b].el==this&&(this.unmousedown(r[b].start),r.splice(b,1),k.unbind("snap.drag.*."+this.id));!r.length&&a.unmousemove(e).unmouseup(l);return this}});C.plugin(function(a,v,y,C){y=y.prototype;var A=/^\s*url\((.+)\)/,w=String,z=a._.$;a.filter={};y.filter=function(d){var f=this;"svg"!=f.type&&(f=f.paper);d=a.parse(w(d));var k=a._.id(),u=z("filter");z(u,{id:k,filterUnits:"userSpaceOnUse"});u.appendChild(d.node);f.defs.appendChild(u);return new v(u)};k.on("snap.util.getattr.filter",
function(){k.stop();var d=z(this.node,"filter");if(d)return(d=w(d).match(A))&&a.select(d[1])});k.on("snap.util.attr.filter",function(d){if(d instanceof v&&"filter"==d.type){k.stop();var f=d.node.id;f||(z(d.node,{id:d.id}),f=d.id);z(this.node,{filter:a.url(f)})}d&&"none"!=d||(k.stop(),this.node.removeAttribute("filter"))});a.filter.blur=function(d,f){null==d&&(d=2);return a.format('<feGaussianBlur stdDeviation="{def}"/>',{def:null==f?d:[d,f]})};a.filter.blur.toString=function(){return this()};a.filter.shadow=
function(d,f,k,u,p){"string"==typeof k&&(p=u=k,k=4);"string"!=typeof u&&(p=u,u="#000");null==k&&(k=4);null==p&&(p=1);null==d&&(d=0,f=2);null==f&&(f=d);u=a.color(u||"#000");return a.format('<feGaussianBlur in="SourceAlpha" stdDeviation="{blur}"/><feOffset dx="{dx}" dy="{dy}" result="offsetblur"/><feFlood flood-color="{color}"/><feComposite in2="offsetblur" operator="in"/><feComponentTransfer><feFuncA type="linear" slope="{opacity}"/></feComponentTransfer><feMerge><feMergeNode/><feMergeNode in="SourceGraphic"/></feMerge>',
{color:u,dx:d,dy:f,blur:k,opacity:p})};a.filter.shadow.toString=function(){return this()};a.filter.grayscale=function(d){null==d&&(d=1);return a.format('<feColorMatrix type="matrix" values="{a} {b} {c} 0 0 {d} {e} {f} 0 0 {g} {b} {h} 0 0 0 0 0 1 0"/>',{a:0.2126+0.7874*(1-d),b:0.7152-0.7152*(1-d),c:0.0722-0.0722*(1-d),d:0.2126-0.2126*(1-d),e:0.7152+0.2848*(1-d),f:0.0722-0.0722*(1-d),g:0.2126-0.2126*(1-d),h:0.0722+0.9278*(1-d)})};a.filter.grayscale.toString=function(){return this()};a.filter.sepia=
function(d){null==d&&(d=1);return a.format('<feColorMatrix type="matrix" values="{a} {b} {c} 0 0 {d} {e} {f} 0 0 {g} {h} {i} 0 0 0 0 0 1 0"/>',{a:0.393+0.607*(1-d),b:0.769-0.769*(1-d),c:0.189-0.189*(1-d),d:0.349-0.349*(1-d),e:0.686+0.314*(1-d),f:0.168-0.168*(1-d),g:0.272-0.272*(1-d),h:0.534-0.534*(1-d),i:0.131+0.869*(1-d)})};a.filter.sepia.toString=function(){return this()};a.filter.saturate=function(d){null==d&&(d=1);return a.format('<feColorMatrix type="saturate" values="{amount}"/>',{amount:1-
d})};a.filter.saturate.toString=function(){return this()};a.filter.hueRotate=function(d){return a.format('<feColorMatrix type="hueRotate" values="{angle}"/>',{angle:d||0})};a.filter.hueRotate.toString=function(){return this()};a.filter.invert=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="table" tableValues="{amount} {amount2}"/><feFuncG type="table" tableValues="{amount} {amount2}"/><feFuncB type="table" tableValues="{amount} {amount2}"/></feComponentTransfer>',{amount:d,
amount2:1-d})};a.filter.invert.toString=function(){return this()};a.filter.brightness=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="linear" slope="{amount}"/><feFuncG type="linear" slope="{amount}"/><feFuncB type="linear" slope="{amount}"/></feComponentTransfer>',{amount:d})};a.filter.brightness.toString=function(){return this()};a.filter.contrast=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="linear" slope="{amount}" intercept="{amount2}"/><feFuncG type="linear" slope="{amount}" intercept="{amount2}"/><feFuncB type="linear" slope="{amount}" intercept="{amount2}"/></feComponentTransfer>',
{amount:d,amount2:0.5-d/2})};a.filter.contrast.toString=function(){return this()}});return C});

]]> </script>
</svg>


<pre class='language-julia'><code class='language-julia'>mat(multi_control)</code></pre>
<pre id='var-hash660344' class='code-output documenter-example-output'>1024×1024 SparseMatrixCSC{ComplexF64, Int64} with 1152 stored entries:
⠑⢄⠀⠀⠀⠀⠀⠲⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠑⢄⠀⠀⠀⠀⠀⠑⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠛⢄⠀⠀⠀⠀⠈⠃⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⢠⡀⠀⠀⠀⠀⠑⢄⠀⠀⠀⠀⠈⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⢄⠀⠀⠀⠀⠀⠑⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠐⠦⠀⠀⠀⠀⠀⠑⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠐⠆⠀⠀⠀⠀⠀⠱⢆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠀⠀⠀⠀⠀⠑⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠀⠀⠀⠀⠈⠃⢀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⡀⠀⠀⠀⠀⠁⣄⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢄⠀⠀⠀⠀⠈⠑⢄⠀⠀⠀⠀⠈⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠦⠀⠀⠀⠀⠀⠱⢆⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⠄⠀⠀⠀⠀⠀⠑⢄⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠂⠀⠀⠀⠀⠀⠑⢄</pre>


<div class="markdown"><p><code>ChainBlock</code> is a composite block for matrix multiplication &#40;reversed order&#41;, it is defined as:</p>
<p class="tex">$$\texttt&#123;chain&#40;G₁, G₂, ..., Gₘ&#41;&#125; &#61; G_m G_&#123;m-1&#125;\ldots G_1
$$</p>
</div>

<pre class='language-julia'><code class='language-julia'>qft3_first = chain(put(3, 1=&gt;H), control(3, 2, 1=&gt;shift(π/2)), control(3, 3, 1=&gt;shift(π/4)))</code></pre>
<pre id='var-qft3_first' class='code-output documenter-example-output'>nqubits: 3
chain
├─ put on (1)
│  └─ H
├─ control(2)
│  └─ (1,) shift(1.5707963267948966)
└─ control(3)
   └─ (1,) shift(0.7853981633974483)
</pre>

<pre class='language-julia'><code class='language-julia'>vizcircuit(qft3_first)</code></pre>
<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:xlink="http://www.w3.org/1999/xlink"
     version="1.2"
     width="85mm" height="56.25mm" viewBox="0 0 85 56.25"
     stroke="none"
     fill="#000000"
     stroke-width="0.3"
     font-size="3.88"

     id="img-06a9486a">
<defs>
  <marker id="arrow" markerWidth="15" markerHeight="7" refX="5" refY="3.5" orient="auto" markerUnits="strokeWidth">
    <path d="M0,0 L15,3.5 L0,7 z" stroke="context-stroke" fill="context-stroke"/>
  </marker>
</defs>
<g stroke-width="0" fill="#000000" id="img-06a9486a-1">
  <g transform="translate(65.62,43.75)">
    <circle cx="0" cy="0" r="0.83" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-06a9486a-2">
  <g transform="translate(65.62,6.25)">
    <path d="M-7.2,-5 L7.2,-5 7.2,5 -7.2,5  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0" fill="#000000" id="img-06a9486a-3">
  <g transform="translate(44.38,25)">
    <circle cx="0" cy="0" r="0.83" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-06a9486a-4">
  <g transform="translate(44.38,6.25)">
    <path d="M-7.2,-5 L7.2,-5 7.2,5 -7.2,5  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-06a9486a-5">
  <g transform="translate(23.12,6.25)">
    <path d="M-5,-5 L5,-5 5,5 -5,5  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-06a9486a-6">
  <g transform="translate(76.25,43.75)">
    <path fill="none" d="M0,0 L0,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-06a9486a-7">
  <g transform="translate(76.25,25)">
    <path fill="none" d="M0,0 L0,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-06a9486a-8">
  <g transform="translate(76.25,6.25)">
    <path fill="none" d="M0,0 L0,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-06a9486a-9">
  <g transform="translate(65.62,25)">
    <path fill="none" d="M-10.62,0 L10.62,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-06a9486a-10">
  <g transform="translate(71.35,43.75)">
    <path fill="none" d="M-4.9,0 L4.9,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-06a9486a-11">
  <g transform="translate(38.65,43.75)">
    <path fill="none" d="M-26.15,0 L26.15,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-06a9486a-12">
  <g transform="translate(74.54,6.25)">
    <path fill="none" d="M-1.71,0 L1.71,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-06a9486a-13">
  <g transform="translate(56.71,6.25)">
    <path fill="none" d="M-1.71,0 L1.71,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-06a9486a-14">
  <g transform="translate(65.62,27.08)">
    <path fill="none" d="M0,15.83 L0,-15.83 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-06a9486a-15">
  <g transform="translate(50.1,25)">
    <path fill="none" d="M-4.9,0 L4.9,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-06a9486a-16">
  <g transform="translate(28.02,25)">
    <path fill="none" d="M-15.52,0 L15.52,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-06a9486a-17">
  <g transform="translate(53.29,6.25)">
    <path fill="none" d="M-1.71,0 L1.71,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-06a9486a-18">
  <g transform="translate(35.46,6.25)">
    <path fill="none" d="M-1.71,0 L1.71,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-06a9486a-19">
  <g transform="translate(44.38,17.71)">
    <path fill="none" d="M0,6.46 L0,-6.46 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-06a9486a-20">
  <g transform="translate(30.94,6.25)">
    <path fill="none" d="M-2.81,0 L2.81,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-06a9486a-21">
  <g transform="translate(15.31,6.25)">
    <path fill="none" d="M-2.81,0 L2.81,0 " class="primitive"/>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="5.64" id="img-06a9486a-22">
  <g transform="translate(65.62,43.75)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em"></text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="3.53" id="img-06a9486a-23">
  <g transform="translate(65.62,6.25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">ϕ(π/4)</text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="5.64" id="img-06a9486a-24">
  <g transform="translate(44.38,25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em"></text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="3.53" id="img-06a9486a-25">
  <g transform="translate(44.38,6.25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">ϕ(π/2)</text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="5.64" id="img-06a9486a-26">
  <g transform="translate(23.12,6.25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">H</text>
    </g>
  </g>
</g>
<script> <![CDATA[
(function(N){var k=/[\.\/]/,L=/\s*,\s*/,C=function(a,d){return a-d},a,v,y={n:{}},M=function(){for(var a=0,d=this.length;a<d;a++)if("undefined"!=typeof this[a])return this[a]},A=function(){for(var a=this.length;--a;)if("undefined"!=typeof this[a])return this[a]},w=function(k,d){k=String(k);var f=v,n=Array.prototype.slice.call(arguments,2),u=w.listeners(k),p=0,b,q=[],e={},l=[],r=a;l.firstDefined=M;l.lastDefined=A;a=k;for(var s=v=0,x=u.length;s<x;s++)"zIndex"in u[s]&&(q.push(u[s].zIndex),0>u[s].zIndex&&
(e[u[s].zIndex]=u[s]));for(q.sort(C);0>q[p];)if(b=e[q[p++] ],l.push(b.apply(d,n)),v)return v=f,l;for(s=0;s<x;s++)if(b=u[s],"zIndex"in b)if(b.zIndex==q[p]){l.push(b.apply(d,n));if(v)break;do if(p++,(b=e[q[p] ])&&l.push(b.apply(d,n)),v)break;while(b)}else e[b.zIndex]=b;else if(l.push(b.apply(d,n)),v)break;v=f;a=r;return l};w._events=y;w.listeners=function(a){a=a.split(k);var d=y,f,n,u,p,b,q,e,l=[d],r=[];u=0;for(p=a.length;u<p;u++){e=[];b=0;for(q=l.length;b<q;b++)for(d=l[b].n,f=[d[a[u] ],d["*"] ],n=2;n--;)if(d=
f[n])e.push(d),r=r.concat(d.f||[]);l=e}return r};w.on=function(a,d){a=String(a);if("function"!=typeof d)return function(){};for(var f=a.split(L),n=0,u=f.length;n<u;n++)(function(a){a=a.split(k);for(var b=y,f,e=0,l=a.length;e<l;e++)b=b.n,b=b.hasOwnProperty(a[e])&&b[a[e] ]||(b[a[e] ]={n:{}});b.f=b.f||[];e=0;for(l=b.f.length;e<l;e++)if(b.f[e]==d){f=!0;break}!f&&b.f.push(d)})(f[n]);return function(a){+a==+a&&(d.zIndex=+a)}};w.f=function(a){var d=[].slice.call(arguments,1);return function(){w.apply(null,
[a,null].concat(d).concat([].slice.call(arguments,0)))}};w.stop=function(){v=1};w.nt=function(k){return k?(new RegExp("(?:\\.|\\/|^)"+k+"(?:\\.|\\/|$)")).test(a):a};w.nts=function(){return a.split(k)};w.off=w.unbind=function(a,d){if(a){var f=a.split(L);if(1<f.length)for(var n=0,u=f.length;n<u;n++)w.off(f[n],d);else{for(var f=a.split(k),p,b,q,e,l=[y],n=0,u=f.length;n<u;n++)for(e=0;e<l.length;e+=q.length-2){q=[e,1];p=l[e].n;if("*"!=f[n])p[f[n] ]&&q.push(p[f[n] ]);else for(b in p)p.hasOwnProperty(b)&&
q.push(p[b]);l.splice.apply(l,q)}n=0;for(u=l.length;n<u;n++)for(p=l[n];p.n;){if(d){if(p.f){e=0;for(f=p.f.length;e<f;e++)if(p.f[e]==d){p.f.splice(e,1);break}!p.f.length&&delete p.f}for(b in p.n)if(p.n.hasOwnProperty(b)&&p.n[b].f){q=p.n[b].f;e=0;for(f=q.length;e<f;e++)if(q[e]==d){q.splice(e,1);break}!q.length&&delete p.n[b].f}}else for(b in delete p.f,p.n)p.n.hasOwnProperty(b)&&p.n[b].f&&delete p.n[b].f;p=p.n}}}else w._events=y={n:{}}};w.once=function(a,d){var f=function(){w.unbind(a,f);return d.apply(this,
arguments)};return w.on(a,f)};w.version="0.4.2";w.toString=function(){return"You are running Eve 0.4.2"};"undefined"!=typeof module&&module.exports?module.exports=w:"function"===typeof define&&define.amd?define("eve",[],function(){return w}):N.eve=w})(this);
(function(N,k){"function"===typeof define&&define.amd?define("Snap.svg",["eve"],function(L){return k(N,L)}):k(N,N.eve)})(this,function(N,k){var L=function(a){var k={},y=N.requestAnimationFrame||N.webkitRequestAnimationFrame||N.mozRequestAnimationFrame||N.oRequestAnimationFrame||N.msRequestAnimationFrame||function(a){setTimeout(a,16)},M=Array.isArray||function(a){return a instanceof Array||"[object Array]"==Object.prototype.toString.call(a)},A=0,w="M"+(+new Date).toString(36),z=function(a){if(null==
a)return this.s;var b=this.s-a;this.b+=this.dur*b;this.B+=this.dur*b;this.s=a},d=function(a){if(null==a)return this.spd;this.spd=a},f=function(a){if(null==a)return this.dur;this.s=this.s*a/this.dur;this.dur=a},n=function(){delete k[this.id];this.update();a("mina.stop."+this.id,this)},u=function(){this.pdif||(delete k[this.id],this.update(),this.pdif=this.get()-this.b)},p=function(){this.pdif&&(this.b=this.get()-this.pdif,delete this.pdif,k[this.id]=this)},b=function(){var a;if(M(this.start)){a=[];
for(var b=0,e=this.start.length;b<e;b++)a[b]=+this.start[b]+(this.end[b]-this.start[b])*this.easing(this.s)}else a=+this.start+(this.end-this.start)*this.easing(this.s);this.set(a)},q=function(){var l=0,b;for(b in k)if(k.hasOwnProperty(b)){var e=k[b],f=e.get();l++;e.s=(f-e.b)/(e.dur/e.spd);1<=e.s&&(delete k[b],e.s=1,l--,function(b){setTimeout(function(){a("mina.finish."+b.id,b)})}(e));e.update()}l&&y(q)},e=function(a,r,s,x,G,h,J){a={id:w+(A++).toString(36),start:a,end:r,b:s,s:0,dur:x-s,spd:1,get:G,
set:h,easing:J||e.linear,status:z,speed:d,duration:f,stop:n,pause:u,resume:p,update:b};k[a.id]=a;r=0;for(var K in k)if(k.hasOwnProperty(K)&&(r++,2==r))break;1==r&&y(q);return a};e.time=Date.now||function(){return+new Date};e.getById=function(a){return k[a]||null};e.linear=function(a){return a};e.easeout=function(a){return Math.pow(a,1.7)};e.easein=function(a){return Math.pow(a,0.48)};e.easeinout=function(a){if(1==a)return 1;if(0==a)return 0;var b=0.48-a/1.04,e=Math.sqrt(0.1734+b*b);a=e-b;a=Math.pow(Math.abs(a),
1/3)*(0>a?-1:1);b=-e-b;b=Math.pow(Math.abs(b),1/3)*(0>b?-1:1);a=a+b+0.5;return 3*(1-a)*a*a+a*a*a};e.backin=function(a){return 1==a?1:a*a*(2.70158*a-1.70158)};e.backout=function(a){if(0==a)return 0;a-=1;return a*a*(2.70158*a+1.70158)+1};e.elastic=function(a){return a==!!a?a:Math.pow(2,-10*a)*Math.sin(2*(a-0.075)*Math.PI/0.3)+1};e.bounce=function(a){a<1/2.75?a*=7.5625*a:a<2/2.75?(a-=1.5/2.75,a=7.5625*a*a+0.75):a<2.5/2.75?(a-=2.25/2.75,a=7.5625*a*a+0.9375):(a-=2.625/2.75,a=7.5625*a*a+0.984375);return a};
return N.mina=e}("undefined"==typeof k?function(){}:k),C=function(){function a(c,t){if(c){if(c.tagName)return x(c);if(y(c,"array")&&a.set)return a.set.apply(a,c);if(c instanceof e)return c;if(null==t)return c=G.doc.querySelector(c),x(c)}return new s(null==c?"100%":c,null==t?"100%":t)}function v(c,a){if(a){"#text"==c&&(c=G.doc.createTextNode(a.text||""));"string"==typeof c&&(c=v(c));if("string"==typeof a)return"xlink:"==a.substring(0,6)?c.getAttributeNS(m,a.substring(6)):"xml:"==a.substring(0,4)?c.getAttributeNS(la,
a.substring(4)):c.getAttribute(a);for(var da in a)if(a[h](da)){var b=J(a[da]);b?"xlink:"==da.substring(0,6)?c.setAttributeNS(m,da.substring(6),b):"xml:"==da.substring(0,4)?c.setAttributeNS(la,da.substring(4),b):c.setAttribute(da,b):c.removeAttribute(da)}}else c=G.doc.createElementNS(la,c);return c}function y(c,a){a=J.prototype.toLowerCase.call(a);return"finite"==a?isFinite(c):"array"==a&&(c instanceof Array||Array.isArray&&Array.isArray(c))?!0:"null"==a&&null===c||a==typeof c&&null!==c||"object"==
a&&c===Object(c)||$.call(c).slice(8,-1).toLowerCase()==a}function M(c){if("function"==typeof c||Object(c)!==c)return c;var a=new c.constructor,b;for(b in c)c[h](b)&&(a[b]=M(c[b]));return a}function A(c,a,b){function m(){var e=Array.prototype.slice.call(arguments,0),f=e.join("\u2400"),d=m.cache=m.cache||{},l=m.count=m.count||[];if(d[h](f)){a:for(var e=l,l=f,B=0,H=e.length;B<H;B++)if(e[B]===l){e.push(e.splice(B,1)[0]);break a}return b?b(d[f]):d[f]}1E3<=l.length&&delete d[l.shift()];l.push(f);d[f]=c.apply(a,
e);return b?b(d[f]):d[f]}return m}function w(c,a,b,m,e,f){return null==e?(c-=b,a-=m,c||a?(180*I.atan2(-a,-c)/C+540)%360:0):w(c,a,e,f)-w(b,m,e,f)}function z(c){return c%360*C/180}function d(c){var a=[];c=c.replace(/(?:^|\s)(\w+)\(([^)]+)\)/g,function(c,b,m){m=m.split(/\s*,\s*|\s+/);"rotate"==b&&1==m.length&&m.push(0,0);"scale"==b&&(2<m.length?m=m.slice(0,2):2==m.length&&m.push(0,0),1==m.length&&m.push(m[0],0,0));"skewX"==b?a.push(["m",1,0,I.tan(z(m[0])),1,0,0]):"skewY"==b?a.push(["m",1,I.tan(z(m[0])),
0,1,0,0]):a.push([b.charAt(0)].concat(m));return c});return a}function f(c,t){var b=O(c),m=new a.Matrix;if(b)for(var e=0,f=b.length;e<f;e++){var h=b[e],d=h.length,B=J(h[0]).toLowerCase(),H=h[0]!=B,l=H?m.invert():0,E;"t"==B&&2==d?m.translate(h[1],0):"t"==B&&3==d?H?(d=l.x(0,0),B=l.y(0,0),H=l.x(h[1],h[2]),l=l.y(h[1],h[2]),m.translate(H-d,l-B)):m.translate(h[1],h[2]):"r"==B?2==d?(E=E||t,m.rotate(h[1],E.x+E.width/2,E.y+E.height/2)):4==d&&(H?(H=l.x(h[2],h[3]),l=l.y(h[2],h[3]),m.rotate(h[1],H,l)):m.rotate(h[1],
h[2],h[3])):"s"==B?2==d||3==d?(E=E||t,m.scale(h[1],h[d-1],E.x+E.width/2,E.y+E.height/2)):4==d?H?(H=l.x(h[2],h[3]),l=l.y(h[2],h[3]),m.scale(h[1],h[1],H,l)):m.scale(h[1],h[1],h[2],h[3]):5==d&&(H?(H=l.x(h[3],h[4]),l=l.y(h[3],h[4]),m.scale(h[1],h[2],H,l)):m.scale(h[1],h[2],h[3],h[4])):"m"==B&&7==d&&m.add(h[1],h[2],h[3],h[4],h[5],h[6])}return m}function n(c,t){if(null==t){var m=!0;t="linearGradient"==c.type||"radialGradient"==c.type?c.node.getAttribute("gradientTransform"):"pattern"==c.type?c.node.getAttribute("patternTransform"):
c.node.getAttribute("transform");if(!t)return new a.Matrix;t=d(t)}else t=a._.rgTransform.test(t)?J(t).replace(/\.{3}|\u2026/g,c._.transform||aa):d(t),y(t,"array")&&(t=a.path?a.path.toString.call(t):J(t)),c._.transform=t;var b=f(t,c.getBBox(1));if(m)return b;c.matrix=b}function u(c){c=c.node.ownerSVGElement&&x(c.node.ownerSVGElement)||c.node.parentNode&&x(c.node.parentNode)||a.select("svg")||a(0,0);var t=c.select("defs"),t=null==t?!1:t.node;t||(t=r("defs",c.node).node);return t}function p(c){return c.node.ownerSVGElement&&
x(c.node.ownerSVGElement)||a.select("svg")}function b(c,a,m){function b(c){if(null==c)return aa;if(c==+c)return c;v(B,{width:c});try{return B.getBBox().width}catch(a){return 0}}function h(c){if(null==c)return aa;if(c==+c)return c;v(B,{height:c});try{return B.getBBox().height}catch(a){return 0}}function e(b,B){null==a?d[b]=B(c.attr(b)||0):b==a&&(d=B(null==m?c.attr(b)||0:m))}var f=p(c).node,d={},B=f.querySelector(".svg---mgr");B||(B=v("rect"),v(B,{x:-9E9,y:-9E9,width:10,height:10,"class":"svg---mgr",
fill:"none"}),f.appendChild(B));switch(c.type){case "rect":e("rx",b),e("ry",h);case "image":e("width",b),e("height",h);case "text":e("x",b);e("y",h);break;case "circle":e("cx",b);e("cy",h);e("r",b);break;case "ellipse":e("cx",b);e("cy",h);e("rx",b);e("ry",h);break;case "line":e("x1",b);e("x2",b);e("y1",h);e("y2",h);break;case "marker":e("refX",b);e("markerWidth",b);e("refY",h);e("markerHeight",h);break;case "radialGradient":e("fx",b);e("fy",h);break;case "tspan":e("dx",b);e("dy",h);break;default:e(a,
b)}f.removeChild(B);return d}function q(c){y(c,"array")||(c=Array.prototype.slice.call(arguments,0));for(var a=0,b=0,m=this.node;this[a];)delete this[a++];for(a=0;a<c.length;a++)"set"==c[a].type?c[a].forEach(function(c){m.appendChild(c.node)}):m.appendChild(c[a].node);for(var h=m.childNodes,a=0;a<h.length;a++)this[b++]=x(h[a]);return this}function e(c){if(c.snap in E)return E[c.snap];var a=this.id=V(),b;try{b=c.ownerSVGElement}catch(m){}this.node=c;b&&(this.paper=new s(b));this.type=c.tagName;this.anims=
{};this._={transform:[]};c.snap=a;E[a]=this;"g"==this.type&&(this.add=q);if(this.type in{g:1,mask:1,pattern:1})for(var e in s.prototype)s.prototype[h](e)&&(this[e]=s.prototype[e])}function l(c){this.node=c}function r(c,a){var b=v(c);a.appendChild(b);return x(b)}function s(c,a){var b,m,f,d=s.prototype;if(c&&"svg"==c.tagName){if(c.snap in E)return E[c.snap];var l=c.ownerDocument;b=new e(c);m=c.getElementsByTagName("desc")[0];f=c.getElementsByTagName("defs")[0];m||(m=v("desc"),m.appendChild(l.createTextNode("Created with Snap")),
b.node.appendChild(m));f||(f=v("defs"),b.node.appendChild(f));b.defs=f;for(var ca in d)d[h](ca)&&(b[ca]=d[ca]);b.paper=b.root=b}else b=r("svg",G.doc.body),v(b.node,{height:a,version:1.1,width:c,xmlns:la});return b}function x(c){return!c||c instanceof e||c instanceof l?c:c.tagName&&"svg"==c.tagName.toLowerCase()?new s(c):c.tagName&&"object"==c.tagName.toLowerCase()&&"image/svg+xml"==c.type?new s(c.contentDocument.getElementsByTagName("svg")[0]):new e(c)}a.version="0.3.0";a.toString=function(){return"Snap v"+
this.version};a._={};var G={win:N,doc:N.document};a._.glob=G;var h="hasOwnProperty",J=String,K=parseFloat,U=parseInt,I=Math,P=I.max,Q=I.min,Y=I.abs,C=I.PI,aa="",$=Object.prototype.toString,F=/^\s*((#[a-f\d]{6})|(#[a-f\d]{3})|rgba?\(\s*([\d\.]+%?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+%?(?:\s*,\s*[\d\.]+%?)?)\s*\)|hsba?\(\s*([\d\.]+(?:deg|\xb0|%)?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+(?:%?\s*,\s*[\d\.]+)?%?)\s*\)|hsla?\(\s*([\d\.]+(?:deg|\xb0|%)?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+(?:%?\s*,\s*[\d\.]+)?%?)\s*\))\s*$/i;a._.separator=
RegExp("[,\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]+");var S=RegExp("[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*"),X={hs:1,rg:1},W=RegExp("([a-z])[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029,]*((-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*)+)",
"ig"),ma=RegExp("([rstm])[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029,]*((-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*)+)","ig"),Z=RegExp("(-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?)[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*",
"ig"),na=0,ba="S"+(+new Date).toString(36),V=function(){return ba+(na++).toString(36)},m="http://www.w3.org/1999/xlink",la="http://www.w3.org/2000/svg",E={},ca=a.url=function(c){return"url('#"+c+"')"};a._.$=v;a._.id=V;a.format=function(){var c=/\{([^\}]+)\}/g,a=/(?:(?:^|\.)(.+?)(?=\[|\.|$|\()|\[('|")(.+?)\2\])(\(\))?/g,b=function(c,b,m){var h=m;b.replace(a,function(c,a,b,m,t){a=a||m;h&&(a in h&&(h=h[a]),"function"==typeof h&&t&&(h=h()))});return h=(null==h||h==m?c:h)+""};return function(a,m){return J(a).replace(c,
function(c,a){return b(c,a,m)})}}();a._.clone=M;a._.cacher=A;a.rad=z;a.deg=function(c){return 180*c/C%360};a.angle=w;a.is=y;a.snapTo=function(c,a,b){b=y(b,"finite")?b:10;if(y(c,"array"))for(var m=c.length;m--;){if(Y(c[m]-a)<=b)return c[m]}else{c=+c;m=a%c;if(m<b)return a-m;if(m>c-b)return a-m+c}return a};a.getRGB=A(function(c){if(!c||(c=J(c)).indexOf("-")+1)return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka};if("none"==c)return{r:-1,g:-1,b:-1,hex:"none",toString:ka};!X[h](c.toLowerCase().substring(0,
2))&&"#"!=c.charAt()&&(c=T(c));if(!c)return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka};var b,m,e,f,d;if(c=c.match(F)){c[2]&&(e=U(c[2].substring(5),16),m=U(c[2].substring(3,5),16),b=U(c[2].substring(1,3),16));c[3]&&(e=U((d=c[3].charAt(3))+d,16),m=U((d=c[3].charAt(2))+d,16),b=U((d=c[3].charAt(1))+d,16));c[4]&&(d=c[4].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b*=2.55),m=K(d[1]),"%"==d[1].slice(-1)&&(m*=2.55),e=K(d[2]),"%"==d[2].slice(-1)&&(e*=2.55),"rgba"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),
d[3]&&"%"==d[3].slice(-1)&&(f/=100));if(c[5])return d=c[5].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b/=100),m=K(d[1]),"%"==d[1].slice(-1)&&(m/=100),e=K(d[2]),"%"==d[2].slice(-1)&&(e/=100),"deg"!=d[0].slice(-3)&&"\u00b0"!=d[0].slice(-1)||(b/=360),"hsba"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),d[3]&&"%"==d[3].slice(-1)&&(f/=100),a.hsb2rgb(b,m,e,f);if(c[6])return d=c[6].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b/=100),m=K(d[1]),"%"==d[1].slice(-1)&&(m/=100),e=K(d[2]),"%"==d[2].slice(-1)&&(e/=100),
"deg"!=d[0].slice(-3)&&"\u00b0"!=d[0].slice(-1)||(b/=360),"hsla"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),d[3]&&"%"==d[3].slice(-1)&&(f/=100),a.hsl2rgb(b,m,e,f);b=Q(I.round(b),255);m=Q(I.round(m),255);e=Q(I.round(e),255);f=Q(P(f,0),1);c={r:b,g:m,b:e,toString:ka};c.hex="#"+(16777216|e|m<<8|b<<16).toString(16).slice(1);c.opacity=y(f,"finite")?f:1;return c}return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka}},a);a.hsb=A(function(c,b,m){return a.hsb2rgb(c,b,m).hex});a.hsl=A(function(c,b,m){return a.hsl2rgb(c,
b,m).hex});a.rgb=A(function(c,a,b,m){if(y(m,"finite")){var e=I.round;return"rgba("+[e(c),e(a),e(b),+m.toFixed(2)]+")"}return"#"+(16777216|b|a<<8|c<<16).toString(16).slice(1)});var T=function(c){var a=G.doc.getElementsByTagName("head")[0]||G.doc.getElementsByTagName("svg")[0];T=A(function(c){if("red"==c.toLowerCase())return"rgb(255, 0, 0)";a.style.color="rgb(255, 0, 0)";a.style.color=c;c=G.doc.defaultView.getComputedStyle(a,aa).getPropertyValue("color");return"rgb(255, 0, 0)"==c?null:c});return T(c)},
qa=function(){return"hsb("+[this.h,this.s,this.b]+")"},ra=function(){return"hsl("+[this.h,this.s,this.l]+")"},ka=function(){return 1==this.opacity||null==this.opacity?this.hex:"rgba("+[this.r,this.g,this.b,this.opacity]+")"},D=function(c,b,m){null==b&&y(c,"object")&&"r"in c&&"g"in c&&"b"in c&&(m=c.b,b=c.g,c=c.r);null==b&&y(c,string)&&(m=a.getRGB(c),c=m.r,b=m.g,m=m.b);if(1<c||1<b||1<m)c/=255,b/=255,m/=255;return[c,b,m]},oa=function(c,b,m,e){c=I.round(255*c);b=I.round(255*b);m=I.round(255*m);c={r:c,
g:b,b:m,opacity:y(e,"finite")?e:1,hex:a.rgb(c,b,m),toString:ka};y(e,"finite")&&(c.opacity=e);return c};a.color=function(c){var b;y(c,"object")&&"h"in c&&"s"in c&&"b"in c?(b=a.hsb2rgb(c),c.r=b.r,c.g=b.g,c.b=b.b,c.opacity=1,c.hex=b.hex):y(c,"object")&&"h"in c&&"s"in c&&"l"in c?(b=a.hsl2rgb(c),c.r=b.r,c.g=b.g,c.b=b.b,c.opacity=1,c.hex=b.hex):(y(c,"string")&&(c=a.getRGB(c)),y(c,"object")&&"r"in c&&"g"in c&&"b"in c&&!("error"in c)?(b=a.rgb2hsl(c),c.h=b.h,c.s=b.s,c.l=b.l,b=a.rgb2hsb(c),c.v=b.b):(c={hex:"none"},
c.r=c.g=c.b=c.h=c.s=c.v=c.l=-1,c.error=1));c.toString=ka;return c};a.hsb2rgb=function(c,a,b,m){y(c,"object")&&"h"in c&&"s"in c&&"b"in c&&(b=c.b,a=c.s,c=c.h,m=c.o);var e,h,d;c=360*c%360/60;d=b*a;a=d*(1-Y(c%2-1));b=e=h=b-d;c=~~c;b+=[d,a,0,0,a,d][c];e+=[a,d,d,a,0,0][c];h+=[0,0,a,d,d,a][c];return oa(b,e,h,m)};a.hsl2rgb=function(c,a,b,m){y(c,"object")&&"h"in c&&"s"in c&&"l"in c&&(b=c.l,a=c.s,c=c.h);if(1<c||1<a||1<b)c/=360,a/=100,b/=100;var e,h,d;c=360*c%360/60;d=2*a*(0.5>b?b:1-b);a=d*(1-Y(c%2-1));b=e=
h=b-d/2;c=~~c;b+=[d,a,0,0,a,d][c];e+=[a,d,d,a,0,0][c];h+=[0,0,a,d,d,a][c];return oa(b,e,h,m)};a.rgb2hsb=function(c,a,b){b=D(c,a,b);c=b[0];a=b[1];b=b[2];var m,e;m=P(c,a,b);e=m-Q(c,a,b);c=((0==e?0:m==c?(a-b)/e:m==a?(b-c)/e+2:(c-a)/e+4)+360)%6*60/360;return{h:c,s:0==e?0:e/m,b:m,toString:qa}};a.rgb2hsl=function(c,a,b){b=D(c,a,b);c=b[0];a=b[1];b=b[2];var m,e,h;m=P(c,a,b);e=Q(c,a,b);h=m-e;c=((0==h?0:m==c?(a-b)/h:m==a?(b-c)/h+2:(c-a)/h+4)+360)%6*60/360;m=(m+e)/2;return{h:c,s:0==h?0:0.5>m?h/(2*m):h/(2-2*
m),l:m,toString:ra}};a.parsePathString=function(c){if(!c)return null;var b=a.path(c);if(b.arr)return a.path.clone(b.arr);var m={a:7,c:6,o:2,h:1,l:2,m:2,r:4,q:4,s:4,t:2,v:1,u:3,z:0},e=[];y(c,"array")&&y(c[0],"array")&&(e=a.path.clone(c));e.length||J(c).replace(W,function(c,a,b){var h=[];c=a.toLowerCase();b.replace(Z,function(c,a){a&&h.push(+a)});"m"==c&&2<h.length&&(e.push([a].concat(h.splice(0,2))),c="l",a="m"==a?"l":"L");"o"==c&&1==h.length&&e.push([a,h[0] ]);if("r"==c)e.push([a].concat(h));else for(;h.length>=
m[c]&&(e.push([a].concat(h.splice(0,m[c]))),m[c]););});e.toString=a.path.toString;b.arr=a.path.clone(e);return e};var O=a.parseTransformString=function(c){if(!c)return null;var b=[];y(c,"array")&&y(c[0],"array")&&(b=a.path.clone(c));b.length||J(c).replace(ma,function(c,a,m){var e=[];a.toLowerCase();m.replace(Z,function(c,a){a&&e.push(+a)});b.push([a].concat(e))});b.toString=a.path.toString;return b};a._.svgTransform2string=d;a._.rgTransform=RegExp("^[a-z][\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*-?\\.?\\d",
"i");a._.transform2matrix=f;a._unit2px=b;a._.getSomeDefs=u;a._.getSomeSVG=p;a.select=function(c){return x(G.doc.querySelector(c))};a.selectAll=function(c){c=G.doc.querySelectorAll(c);for(var b=(a.set||Array)(),m=0;m<c.length;m++)b.push(x(c[m]));return b};setInterval(function(){for(var c in E)if(E[h](c)){var a=E[c],b=a.node;("svg"!=a.type&&!b.ownerSVGElement||"svg"==a.type&&(!b.parentNode||"ownerSVGElement"in b.parentNode&&!b.ownerSVGElement))&&delete E[c]}},1E4);(function(c){function m(c){function a(c,
b){var m=v(c.node,b);(m=(m=m&&m.match(d))&&m[2])&&"#"==m.charAt()&&(m=m.substring(1))&&(f[m]=(f[m]||[]).concat(function(a){var m={};m[b]=ca(a);v(c.node,m)}))}function b(c){var a=v(c.node,"xlink:href");a&&"#"==a.charAt()&&(a=a.substring(1))&&(f[a]=(f[a]||[]).concat(function(a){c.attr("xlink:href","#"+a)}))}var e=c.selectAll("*"),h,d=/^\s*url\(("|'|)(.*)\1\)\s*$/;c=[];for(var f={},l=0,E=e.length;l<E;l++){h=e[l];a(h,"fill");a(h,"stroke");a(h,"filter");a(h,"mask");a(h,"clip-path");b(h);var t=v(h.node,
"id");t&&(v(h.node,{id:h.id}),c.push({old:t,id:h.id}))}l=0;for(E=c.length;l<E;l++)if(e=f[c[l].old])for(h=0,t=e.length;h<t;h++)e[h](c[l].id)}function e(c,a,b){return function(m){m=m.slice(c,a);1==m.length&&(m=m[0]);return b?b(m):m}}function d(c){return function(){var a=c?"<"+this.type:"",b=this.node.attributes,m=this.node.childNodes;if(c)for(var e=0,h=b.length;e<h;e++)a+=" "+b[e].name+'="'+b[e].value.replace(/"/g,'\\"')+'"';if(m.length){c&&(a+=">");e=0;for(h=m.length;e<h;e++)3==m[e].nodeType?a+=m[e].nodeValue:
1==m[e].nodeType&&(a+=x(m[e]).toString());c&&(a+="</"+this.type+">")}else c&&(a+="/>");return a}}c.attr=function(c,a){if(!c)return this;if(y(c,"string"))if(1<arguments.length){var b={};b[c]=a;c=b}else return k("snap.util.getattr."+c,this).firstDefined();for(var m in c)c[h](m)&&k("snap.util.attr."+m,this,c[m]);return this};c.getBBox=function(c){if(!a.Matrix||!a.path)return this.node.getBBox();var b=this,m=new a.Matrix;if(b.removed)return a._.box();for(;"use"==b.type;)if(c||(m=m.add(b.transform().localMatrix.translate(b.attr("x")||
0,b.attr("y")||0))),b.original)b=b.original;else var e=b.attr("xlink:href"),b=b.original=b.node.ownerDocument.getElementById(e.substring(e.indexOf("#")+1));var e=b._,h=a.path.get[b.type]||a.path.get.deflt;try{if(c)return e.bboxwt=h?a.path.getBBox(b.realPath=h(b)):a._.box(b.node.getBBox()),a._.box(e.bboxwt);b.realPath=h(b);b.matrix=b.transform().localMatrix;e.bbox=a.path.getBBox(a.path.map(b.realPath,m.add(b.matrix)));return a._.box(e.bbox)}catch(d){return a._.box()}};var f=function(){return this.string};
c.transform=function(c){var b=this._;if(null==c){var m=this;c=new a.Matrix(this.node.getCTM());for(var e=n(this),h=[e],d=new a.Matrix,l=e.toTransformString(),b=J(e)==J(this.matrix)?J(b.transform):l;"svg"!=m.type&&(m=m.parent());)h.push(n(m));for(m=h.length;m--;)d.add(h[m]);return{string:b,globalMatrix:c,totalMatrix:d,localMatrix:e,diffMatrix:c.clone().add(e.invert()),global:c.toTransformString(),total:d.toTransformString(),local:l,toString:f}}c instanceof a.Matrix?this.matrix=c:n(this,c);this.node&&
("linearGradient"==this.type||"radialGradient"==this.type?v(this.node,{gradientTransform:this.matrix}):"pattern"==this.type?v(this.node,{patternTransform:this.matrix}):v(this.node,{transform:this.matrix}));return this};c.parent=function(){return x(this.node.parentNode)};c.append=c.add=function(c){if(c){if("set"==c.type){var a=this;c.forEach(function(c){a.add(c)});return this}c=x(c);this.node.appendChild(c.node);c.paper=this.paper}return this};c.appendTo=function(c){c&&(c=x(c),c.append(this));return this};
c.prepend=function(c){if(c){if("set"==c.type){var a=this,b;c.forEach(function(c){b?b.after(c):a.prepend(c);b=c});return this}c=x(c);var m=c.parent();this.node.insertBefore(c.node,this.node.firstChild);this.add&&this.add();c.paper=this.paper;this.parent()&&this.parent().add();m&&m.add()}return this};c.prependTo=function(c){c=x(c);c.prepend(this);return this};c.before=function(c){if("set"==c.type){var a=this;c.forEach(function(c){var b=c.parent();a.node.parentNode.insertBefore(c.node,a.node);b&&b.add()});
this.parent().add();return this}c=x(c);var b=c.parent();this.node.parentNode.insertBefore(c.node,this.node);this.parent()&&this.parent().add();b&&b.add();c.paper=this.paper;return this};c.after=function(c){c=x(c);var a=c.parent();this.node.nextSibling?this.node.parentNode.insertBefore(c.node,this.node.nextSibling):this.node.parentNode.appendChild(c.node);this.parent()&&this.parent().add();a&&a.add();c.paper=this.paper;return this};c.insertBefore=function(c){c=x(c);var a=this.parent();c.node.parentNode.insertBefore(this.node,
c.node);this.paper=c.paper;a&&a.add();c.parent()&&c.parent().add();return this};c.insertAfter=function(c){c=x(c);var a=this.parent();c.node.parentNode.insertBefore(this.node,c.node.nextSibling);this.paper=c.paper;a&&a.add();c.parent()&&c.parent().add();return this};c.remove=function(){var c=this.parent();this.node.parentNode&&this.node.parentNode.removeChild(this.node);delete this.paper;this.removed=!0;c&&c.add();return this};c.select=function(c){return x(this.node.querySelector(c))};c.selectAll=
function(c){c=this.node.querySelectorAll(c);for(var b=(a.set||Array)(),m=0;m<c.length;m++)b.push(x(c[m]));return b};c.asPX=function(c,a){null==a&&(a=this.attr(c));return+b(this,c,a)};c.use=function(){var c,a=this.node.id;a||(a=this.id,v(this.node,{id:a}));c="linearGradient"==this.type||"radialGradient"==this.type||"pattern"==this.type?r(this.type,this.node.parentNode):r("use",this.node.parentNode);v(c.node,{"xlink:href":"#"+a});c.original=this;return c};var l=/\S+/g;c.addClass=function(c){var a=(c||
"").match(l)||[];c=this.node;var b=c.className.baseVal,m=b.match(l)||[],e,h,d;if(a.length){for(e=0;d=a[e++];)h=m.indexOf(d),~h||m.push(d);a=m.join(" ");b!=a&&(c.className.baseVal=a)}return this};c.removeClass=function(c){var a=(c||"").match(l)||[];c=this.node;var b=c.className.baseVal,m=b.match(l)||[],e,h;if(m.length){for(e=0;h=a[e++];)h=m.indexOf(h),~h&&m.splice(h,1);a=m.join(" ");b!=a&&(c.className.baseVal=a)}return this};c.hasClass=function(c){return!!~(this.node.className.baseVal.match(l)||[]).indexOf(c)};
c.toggleClass=function(c,a){if(null!=a)return a?this.addClass(c):this.removeClass(c);var b=(c||"").match(l)||[],m=this.node,e=m.className.baseVal,h=e.match(l)||[],d,f,E;for(d=0;E=b[d++];)f=h.indexOf(E),~f?h.splice(f,1):h.push(E);b=h.join(" ");e!=b&&(m.className.baseVal=b);return this};c.clone=function(){var c=x(this.node.cloneNode(!0));v(c.node,"id")&&v(c.node,{id:c.id});m(c);c.insertAfter(this);return c};c.toDefs=function(){u(this).appendChild(this.node);return this};c.pattern=c.toPattern=function(c,
a,b,m){var e=r("pattern",u(this));null==c&&(c=this.getBBox());y(c,"object")&&"x"in c&&(a=c.y,b=c.width,m=c.height,c=c.x);v(e.node,{x:c,y:a,width:b,height:m,patternUnits:"userSpaceOnUse",id:e.id,viewBox:[c,a,b,m].join(" ")});e.node.appendChild(this.node);return e};c.marker=function(c,a,b,m,e,h){var d=r("marker",u(this));null==c&&(c=this.getBBox());y(c,"object")&&"x"in c&&(a=c.y,b=c.width,m=c.height,e=c.refX||c.cx,h=c.refY||c.cy,c=c.x);v(d.node,{viewBox:[c,a,b,m].join(" "),markerWidth:b,markerHeight:m,
orient:"auto",refX:e||0,refY:h||0,id:d.id});d.node.appendChild(this.node);return d};var E=function(c,a,b,m){"function"!=typeof b||b.length||(m=b,b=L.linear);this.attr=c;this.dur=a;b&&(this.easing=b);m&&(this.callback=m)};a._.Animation=E;a.animation=function(c,a,b,m){return new E(c,a,b,m)};c.inAnim=function(){var c=[],a;for(a in this.anims)this.anims[h](a)&&function(a){c.push({anim:new E(a._attrs,a.dur,a.easing,a._callback),mina:a,curStatus:a.status(),status:function(c){return a.status(c)},stop:function(){a.stop()}})}(this.anims[a]);
return c};a.animate=function(c,a,b,m,e,h){"function"!=typeof e||e.length||(h=e,e=L.linear);var d=L.time();c=L(c,a,d,d+m,L.time,b,e);h&&k.once("mina.finish."+c.id,h);return c};c.stop=function(){for(var c=this.inAnim(),a=0,b=c.length;a<b;a++)c[a].stop();return this};c.animate=function(c,a,b,m){"function"!=typeof b||b.length||(m=b,b=L.linear);c instanceof E&&(m=c.callback,b=c.easing,a=b.dur,c=c.attr);var d=[],f=[],l={},t,ca,n,T=this,q;for(q in c)if(c[h](q)){T.equal?(n=T.equal(q,J(c[q])),t=n.from,ca=
n.to,n=n.f):(t=+T.attr(q),ca=+c[q]);var la=y(t,"array")?t.length:1;l[q]=e(d.length,d.length+la,n);d=d.concat(t);f=f.concat(ca)}t=L.time();var p=L(d,f,t,t+a,L.time,function(c){var a={},b;for(b in l)l[h](b)&&(a[b]=l[b](c));T.attr(a)},b);T.anims[p.id]=p;p._attrs=c;p._callback=m;k("snap.animcreated."+T.id,p);k.once("mina.finish."+p.id,function(){delete T.anims[p.id];m&&m.call(T)});k.once("mina.stop."+p.id,function(){delete T.anims[p.id]});return T};var T={};c.data=function(c,b){var m=T[this.id]=T[this.id]||
{};if(0==arguments.length)return k("snap.data.get."+this.id,this,m,null),m;if(1==arguments.length){if(a.is(c,"object")){for(var e in c)c[h](e)&&this.data(e,c[e]);return this}k("snap.data.get."+this.id,this,m[c],c);return m[c]}m[c]=b;k("snap.data.set."+this.id,this,b,c);return this};c.removeData=function(c){null==c?T[this.id]={}:T[this.id]&&delete T[this.id][c];return this};c.outerSVG=c.toString=d(1);c.innerSVG=d()})(e.prototype);a.parse=function(c){var a=G.doc.createDocumentFragment(),b=!0,m=G.doc.createElement("div");
c=J(c);c.match(/^\s*<\s*svg(?:\s|>)/)||(c="<svg>"+c+"</svg>",b=!1);m.innerHTML=c;if(c=m.getElementsByTagName("svg")[0])if(b)a=c;else for(;c.firstChild;)a.appendChild(c.firstChild);m.innerHTML=aa;return new l(a)};l.prototype.select=e.prototype.select;l.prototype.selectAll=e.prototype.selectAll;a.fragment=function(){for(var c=Array.prototype.slice.call(arguments,0),b=G.doc.createDocumentFragment(),m=0,e=c.length;m<e;m++){var h=c[m];h.node&&h.node.nodeType&&b.appendChild(h.node);h.nodeType&&b.appendChild(h);
"string"==typeof h&&b.appendChild(a.parse(h).node)}return new l(b)};a._.make=r;a._.wrap=x;s.prototype.el=function(c,a){var b=r(c,this.node);a&&b.attr(a);return b};k.on("snap.util.getattr",function(){var c=k.nt(),c=c.substring(c.lastIndexOf(".")+1),a=c.replace(/[A-Z]/g,function(c){return"-"+c.toLowerCase()});return pa[h](a)?this.node.ownerDocument.defaultView.getComputedStyle(this.node,null).getPropertyValue(a):v(this.node,c)});var pa={"alignment-baseline":0,"baseline-shift":0,clip:0,"clip-path":0,
"clip-rule":0,color:0,"color-interpolation":0,"color-interpolation-filters":0,"color-profile":0,"color-rendering":0,cursor:0,direction:0,display:0,"dominant-baseline":0,"enable-background":0,fill:0,"fill-opacity":0,"fill-rule":0,filter:0,"flood-color":0,"flood-opacity":0,font:0,"font-family":0,"font-size":0,"font-size-adjust":0,"font-stretch":0,"font-style":0,"font-variant":0,"font-weight":0,"glyph-orientation-horizontal":0,"glyph-orientation-vertical":0,"image-rendering":0,kerning:0,"letter-spacing":0,
"lighting-color":0,marker:0,"marker-end":0,"marker-mid":0,"marker-start":0,mask:0,opacity:0,overflow:0,"pointer-events":0,"shape-rendering":0,"stop-color":0,"stop-opacity":0,stroke:0,"stroke-dasharray":0,"stroke-dashoffset":0,"stroke-linecap":0,"stroke-linejoin":0,"stroke-miterlimit":0,"stroke-opacity":0,"stroke-width":0,"text-anchor":0,"text-decoration":0,"text-rendering":0,"unicode-bidi":0,visibility:0,"word-spacing":0,"writing-mode":0};k.on("snap.util.attr",function(c){var a=k.nt(),b={},a=a.substring(a.lastIndexOf(".")+
1);b[a]=c;var m=a.replace(/-(\w)/gi,function(c,a){return a.toUpperCase()}),a=a.replace(/[A-Z]/g,function(c){return"-"+c.toLowerCase()});pa[h](a)?this.node.style[m]=null==c?aa:c:v(this.node,b)});a.ajax=function(c,a,b,m){var e=new XMLHttpRequest,h=V();if(e){if(y(a,"function"))m=b,b=a,a=null;else if(y(a,"object")){var d=[],f;for(f in a)a.hasOwnProperty(f)&&d.push(encodeURIComponent(f)+"="+encodeURIComponent(a[f]));a=d.join("&")}e.open(a?"POST":"GET",c,!0);a&&(e.setRequestHeader("X-Requested-With","XMLHttpRequest"),
e.setRequestHeader("Content-type","application/x-www-form-urlencoded"));b&&(k.once("snap.ajax."+h+".0",b),k.once("snap.ajax."+h+".200",b),k.once("snap.ajax."+h+".304",b));e.onreadystatechange=function(){4==e.readyState&&k("snap.ajax."+h+"."+e.status,m,e)};if(4==e.readyState)return e;e.send(a);return e}};a.load=function(c,b,m){a.ajax(c,function(c){c=a.parse(c.responseText);m?b.call(m,c):b(c)})};a.getElementByPoint=function(c,a){var b,m,e=G.doc.elementFromPoint(c,a);if(G.win.opera&&"svg"==e.tagName){b=
e;m=b.getBoundingClientRect();b=b.ownerDocument;var h=b.body,d=b.documentElement;b=m.top+(g.win.pageYOffset||d.scrollTop||h.scrollTop)-(d.clientTop||h.clientTop||0);m=m.left+(g.win.pageXOffset||d.scrollLeft||h.scrollLeft)-(d.clientLeft||h.clientLeft||0);h=e.createSVGRect();h.x=c-m;h.y=a-b;h.width=h.height=1;b=e.getIntersectionList(h,null);b.length&&(e=b[b.length-1])}return e?x(e):null};a.plugin=function(c){c(a,e,s,G,l)};return G.win.Snap=a}();C.plugin(function(a,k,y,M,A){function w(a,d,f,b,q,e){null==
d&&"[object SVGMatrix]"==z.call(a)?(this.a=a.a,this.b=a.b,this.c=a.c,this.d=a.d,this.e=a.e,this.f=a.f):null!=a?(this.a=+a,this.b=+d,this.c=+f,this.d=+b,this.e=+q,this.f=+e):(this.a=1,this.c=this.b=0,this.d=1,this.f=this.e=0)}var z=Object.prototype.toString,d=String,f=Math;(function(n){function k(a){return a[0]*a[0]+a[1]*a[1]}function p(a){var d=f.sqrt(k(a));a[0]&&(a[0]/=d);a[1]&&(a[1]/=d)}n.add=function(a,d,e,f,n,p){var k=[[],[],[] ],u=[[this.a,this.c,this.e],[this.b,this.d,this.f],[0,0,1] ];d=[[a,
e,n],[d,f,p],[0,0,1] ];a&&a instanceof w&&(d=[[a.a,a.c,a.e],[a.b,a.d,a.f],[0,0,1] ]);for(a=0;3>a;a++)for(e=0;3>e;e++){for(f=n=0;3>f;f++)n+=u[a][f]*d[f][e];k[a][e]=n}this.a=k[0][0];this.b=k[1][0];this.c=k[0][1];this.d=k[1][1];this.e=k[0][2];this.f=k[1][2];return this};n.invert=function(){var a=this.a*this.d-this.b*this.c;return new w(this.d/a,-this.b/a,-this.c/a,this.a/a,(this.c*this.f-this.d*this.e)/a,(this.b*this.e-this.a*this.f)/a)};n.clone=function(){return new w(this.a,this.b,this.c,this.d,this.e,
this.f)};n.translate=function(a,d){return this.add(1,0,0,1,a,d)};n.scale=function(a,d,e,f){null==d&&(d=a);(e||f)&&this.add(1,0,0,1,e,f);this.add(a,0,0,d,0,0);(e||f)&&this.add(1,0,0,1,-e,-f);return this};n.rotate=function(b,d,e){b=a.rad(b);d=d||0;e=e||0;var l=+f.cos(b).toFixed(9);b=+f.sin(b).toFixed(9);this.add(l,b,-b,l,d,e);return this.add(1,0,0,1,-d,-e)};n.x=function(a,d){return a*this.a+d*this.c+this.e};n.y=function(a,d){return a*this.b+d*this.d+this.f};n.get=function(a){return+this[d.fromCharCode(97+
a)].toFixed(4)};n.toString=function(){return"matrix("+[this.get(0),this.get(1),this.get(2),this.get(3),this.get(4),this.get(5)].join()+")"};n.offset=function(){return[this.e.toFixed(4),this.f.toFixed(4)]};n.determinant=function(){return this.a*this.d-this.b*this.c};n.split=function(){var b={};b.dx=this.e;b.dy=this.f;var d=[[this.a,this.c],[this.b,this.d] ];b.scalex=f.sqrt(k(d[0]));p(d[0]);b.shear=d[0][0]*d[1][0]+d[0][1]*d[1][1];d[1]=[d[1][0]-d[0][0]*b.shear,d[1][1]-d[0][1]*b.shear];b.scaley=f.sqrt(k(d[1]));
p(d[1]);b.shear/=b.scaley;0>this.determinant()&&(b.scalex=-b.scalex);var e=-d[0][1],d=d[1][1];0>d?(b.rotate=a.deg(f.acos(d)),0>e&&(b.rotate=360-b.rotate)):b.rotate=a.deg(f.asin(e));b.isSimple=!+b.shear.toFixed(9)&&(b.scalex.toFixed(9)==b.scaley.toFixed(9)||!b.rotate);b.isSuperSimple=!+b.shear.toFixed(9)&&b.scalex.toFixed(9)==b.scaley.toFixed(9)&&!b.rotate;b.noRotation=!+b.shear.toFixed(9)&&!b.rotate;return b};n.toTransformString=function(a){a=a||this.split();if(+a.shear.toFixed(9))return"m"+[this.get(0),
this.get(1),this.get(2),this.get(3),this.get(4),this.get(5)];a.scalex=+a.scalex.toFixed(4);a.scaley=+a.scaley.toFixed(4);a.rotate=+a.rotate.toFixed(4);return(a.dx||a.dy?"t"+[+a.dx.toFixed(4),+a.dy.toFixed(4)]:"")+(1!=a.scalex||1!=a.scaley?"s"+[a.scalex,a.scaley,0,0]:"")+(a.rotate?"r"+[+a.rotate.toFixed(4),0,0]:"")}})(w.prototype);a.Matrix=w;a.matrix=function(a,d,f,b,k,e){return new w(a,d,f,b,k,e)}});C.plugin(function(a,v,y,M,A){function w(h){return function(d){k.stop();d instanceof A&&1==d.node.childNodes.length&&
("radialGradient"==d.node.firstChild.tagName||"linearGradient"==d.node.firstChild.tagName||"pattern"==d.node.firstChild.tagName)&&(d=d.node.firstChild,b(this).appendChild(d),d=u(d));if(d instanceof v)if("radialGradient"==d.type||"linearGradient"==d.type||"pattern"==d.type){d.node.id||e(d.node,{id:d.id});var f=l(d.node.id)}else f=d.attr(h);else f=a.color(d),f.error?(f=a(b(this).ownerSVGElement).gradient(d))?(f.node.id||e(f.node,{id:f.id}),f=l(f.node.id)):f=d:f=r(f);d={};d[h]=f;e(this.node,d);this.node.style[h]=
x}}function z(a){k.stop();a==+a&&(a+="px");this.node.style.fontSize=a}function d(a){var b=[];a=a.childNodes;for(var e=0,f=a.length;e<f;e++){var l=a[e];3==l.nodeType&&b.push(l.nodeValue);"tspan"==l.tagName&&(1==l.childNodes.length&&3==l.firstChild.nodeType?b.push(l.firstChild.nodeValue):b.push(d(l)))}return b}function f(){k.stop();return this.node.style.fontSize}var n=a._.make,u=a._.wrap,p=a.is,b=a._.getSomeDefs,q=/^url\(#?([^)]+)\)$/,e=a._.$,l=a.url,r=String,s=a._.separator,x="";k.on("snap.util.attr.mask",
function(a){if(a instanceof v||a instanceof A){k.stop();a instanceof A&&1==a.node.childNodes.length&&(a=a.node.firstChild,b(this).appendChild(a),a=u(a));if("mask"==a.type)var d=a;else d=n("mask",b(this)),d.node.appendChild(a.node);!d.node.id&&e(d.node,{id:d.id});e(this.node,{mask:l(d.id)})}});(function(a){k.on("snap.util.attr.clip",a);k.on("snap.util.attr.clip-path",a);k.on("snap.util.attr.clipPath",a)})(function(a){if(a instanceof v||a instanceof A){k.stop();if("clipPath"==a.type)var d=a;else d=
n("clipPath",b(this)),d.node.appendChild(a.node),!d.node.id&&e(d.node,{id:d.id});e(this.node,{"clip-path":l(d.id)})}});k.on("snap.util.attr.fill",w("fill"));k.on("snap.util.attr.stroke",w("stroke"));var G=/^([lr])(?:\(([^)]*)\))?(.*)$/i;k.on("snap.util.grad.parse",function(a){a=r(a);var b=a.match(G);if(!b)return null;a=b[1];var e=b[2],b=b[3],e=e.split(/\s*,\s*/).map(function(a){return+a==a?+a:a});1==e.length&&0==e[0]&&(e=[]);b=b.split("-");b=b.map(function(a){a=a.split(":");var b={color:a[0]};a[1]&&
(b.offset=parseFloat(a[1]));return b});return{type:a,params:e,stops:b}});k.on("snap.util.attr.d",function(b){k.stop();p(b,"array")&&p(b[0],"array")&&(b=a.path.toString.call(b));b=r(b);b.match(/[ruo]/i)&&(b=a.path.toAbsolute(b));e(this.node,{d:b})})(-1);k.on("snap.util.attr.#text",function(a){k.stop();a=r(a);for(a=M.doc.createTextNode(a);this.node.firstChild;)this.node.removeChild(this.node.firstChild);this.node.appendChild(a)})(-1);k.on("snap.util.attr.path",function(a){k.stop();this.attr({d:a})})(-1);
k.on("snap.util.attr.class",function(a){k.stop();this.node.className.baseVal=a})(-1);k.on("snap.util.attr.viewBox",function(a){a=p(a,"object")&&"x"in a?[a.x,a.y,a.width,a.height].join(" "):p(a,"array")?a.join(" "):a;e(this.node,{viewBox:a});k.stop()})(-1);k.on("snap.util.attr.transform",function(a){this.transform(a);k.stop()})(-1);k.on("snap.util.attr.r",function(a){"rect"==this.type&&(k.stop(),e(this.node,{rx:a,ry:a}))})(-1);k.on("snap.util.attr.textpath",function(a){k.stop();if("text"==this.type){var d,
f;if(!a&&this.textPath){for(a=this.textPath;a.node.firstChild;)this.node.appendChild(a.node.firstChild);a.remove();delete this.textPath}else if(p(a,"string")?(d=b(this),a=u(d.parentNode).path(a),d.appendChild(a.node),d=a.id,a.attr({id:d})):(a=u(a),a instanceof v&&(d=a.attr("id"),d||(d=a.id,a.attr({id:d})))),d)if(a=this.textPath,f=this.node,a)a.attr({"xlink:href":"#"+d});else{for(a=e("textPath",{"xlink:href":"#"+d});f.firstChild;)a.appendChild(f.firstChild);f.appendChild(a);this.textPath=u(a)}}})(-1);
k.on("snap.util.attr.text",function(a){if("text"==this.type){for(var b=this.node,d=function(a){var b=e("tspan");if(p(a,"array"))for(var f=0;f<a.length;f++)b.appendChild(d(a[f]));else b.appendChild(M.doc.createTextNode(a));b.normalize&&b.normalize();return b};b.firstChild;)b.removeChild(b.firstChild);for(a=d(a);a.firstChild;)b.appendChild(a.firstChild)}k.stop()})(-1);k.on("snap.util.attr.fontSize",z)(-1);k.on("snap.util.attr.font-size",z)(-1);k.on("snap.util.getattr.transform",function(){k.stop();
return this.transform()})(-1);k.on("snap.util.getattr.textpath",function(){k.stop();return this.textPath})(-1);(function(){function b(d){return function(){k.stop();var b=M.doc.defaultView.getComputedStyle(this.node,null).getPropertyValue("marker-"+d);return"none"==b?b:a(M.doc.getElementById(b.match(q)[1]))}}function d(a){return function(b){k.stop();var d="marker"+a.charAt(0).toUpperCase()+a.substring(1);if(""==b||!b)this.node.style[d]="none";else if("marker"==b.type){var f=b.node.id;f||e(b.node,{id:b.id});
this.node.style[d]=l(f)}}}k.on("snap.util.getattr.marker-end",b("end"))(-1);k.on("snap.util.getattr.markerEnd",b("end"))(-1);k.on("snap.util.getattr.marker-start",b("start"))(-1);k.on("snap.util.getattr.markerStart",b("start"))(-1);k.on("snap.util.getattr.marker-mid",b("mid"))(-1);k.on("snap.util.getattr.markerMid",b("mid"))(-1);k.on("snap.util.attr.marker-end",d("end"))(-1);k.on("snap.util.attr.markerEnd",d("end"))(-1);k.on("snap.util.attr.marker-start",d("start"))(-1);k.on("snap.util.attr.markerStart",
d("start"))(-1);k.on("snap.util.attr.marker-mid",d("mid"))(-1);k.on("snap.util.attr.markerMid",d("mid"))(-1)})();k.on("snap.util.getattr.r",function(){if("rect"==this.type&&e(this.node,"rx")==e(this.node,"ry"))return k.stop(),e(this.node,"rx")})(-1);k.on("snap.util.getattr.text",function(){if("text"==this.type||"tspan"==this.type){k.stop();var a=d(this.node);return 1==a.length?a[0]:a}})(-1);k.on("snap.util.getattr.#text",function(){return this.node.textContent})(-1);k.on("snap.util.getattr.viewBox",
function(){k.stop();var b=e(this.node,"viewBox");if(b)return b=b.split(s),a._.box(+b[0],+b[1],+b[2],+b[3])})(-1);k.on("snap.util.getattr.points",function(){var a=e(this.node,"points");k.stop();if(a)return a.split(s)})(-1);k.on("snap.util.getattr.path",function(){var a=e(this.node,"d");k.stop();return a})(-1);k.on("snap.util.getattr.class",function(){return this.node.className.baseVal})(-1);k.on("snap.util.getattr.fontSize",f)(-1);k.on("snap.util.getattr.font-size",f)(-1)});C.plugin(function(a,v,y,
M,A){function w(a){return a}function z(a){return function(b){return+b.toFixed(3)+a}}var d={"+":function(a,b){return a+b},"-":function(a,b){return a-b},"/":function(a,b){return a/b},"*":function(a,b){return a*b}},f=String,n=/[a-z]+$/i,u=/^\s*([+\-\/*])\s*=\s*([\d.eE+\-]+)\s*([^\d\s]+)?\s*$/;k.on("snap.util.attr",function(a){if(a=f(a).match(u)){var b=k.nt(),b=b.substring(b.lastIndexOf(".")+1),q=this.attr(b),e={};k.stop();var l=a[3]||"",r=q.match(n),s=d[a[1] ];r&&r==l?a=s(parseFloat(q),+a[2]):(q=this.asPX(b),
a=s(this.asPX(b),this.asPX(b,a[2]+l)));isNaN(q)||isNaN(a)||(e[b]=a,this.attr(e))}})(-10);k.on("snap.util.equal",function(a,b){var q=f(this.attr(a)||""),e=f(b).match(u);if(e){k.stop();var l=e[3]||"",r=q.match(n),s=d[e[1] ];if(r&&r==l)return{from:parseFloat(q),to:s(parseFloat(q),+e[2]),f:z(r)};q=this.asPX(a);return{from:q,to:s(q,this.asPX(a,e[2]+l)),f:w}}})(-10)});C.plugin(function(a,v,y,M,A){var w=y.prototype,z=a.is;w.rect=function(a,d,k,p,b,q){var e;null==q&&(q=b);z(a,"object")&&"[object Object]"==
a?e=a:null!=a&&(e={x:a,y:d,width:k,height:p},null!=b&&(e.rx=b,e.ry=q));return this.el("rect",e)};w.circle=function(a,d,k){var p;z(a,"object")&&"[object Object]"==a?p=a:null!=a&&(p={cx:a,cy:d,r:k});return this.el("circle",p)};var d=function(){function a(){this.parentNode.removeChild(this)}return function(d,k){var p=M.doc.createElement("img"),b=M.doc.body;p.style.cssText="position:absolute;left:-9999em;top:-9999em";p.onload=function(){k.call(p);p.onload=p.onerror=null;b.removeChild(p)};p.onerror=a;
b.appendChild(p);p.src=d}}();w.image=function(f,n,k,p,b){var q=this.el("image");if(z(f,"object")&&"src"in f)q.attr(f);else if(null!=f){var e={"xlink:href":f,preserveAspectRatio:"none"};null!=n&&null!=k&&(e.x=n,e.y=k);null!=p&&null!=b?(e.width=p,e.height=b):d(f,function(){a._.$(q.node,{width:this.offsetWidth,height:this.offsetHeight})});a._.$(q.node,e)}return q};w.ellipse=function(a,d,k,p){var b;z(a,"object")&&"[object Object]"==a?b=a:null!=a&&(b={cx:a,cy:d,rx:k,ry:p});return this.el("ellipse",b)};
w.path=function(a){var d;z(a,"object")&&!z(a,"array")?d=a:a&&(d={d:a});return this.el("path",d)};w.group=w.g=function(a){var d=this.el("g");1==arguments.length&&a&&!a.type?d.attr(a):arguments.length&&d.add(Array.prototype.slice.call(arguments,0));return d};w.svg=function(a,d,k,p,b,q,e,l){var r={};z(a,"object")&&null==d?r=a:(null!=a&&(r.x=a),null!=d&&(r.y=d),null!=k&&(r.width=k),null!=p&&(r.height=p),null!=b&&null!=q&&null!=e&&null!=l&&(r.viewBox=[b,q,e,l]));return this.el("svg",r)};w.mask=function(a){var d=
this.el("mask");1==arguments.length&&a&&!a.type?d.attr(a):arguments.length&&d.add(Array.prototype.slice.call(arguments,0));return d};w.ptrn=function(a,d,k,p,b,q,e,l){if(z(a,"object"))var r=a;else arguments.length?(r={},null!=a&&(r.x=a),null!=d&&(r.y=d),null!=k&&(r.width=k),null!=p&&(r.height=p),null!=b&&null!=q&&null!=e&&null!=l&&(r.viewBox=[b,q,e,l])):r={patternUnits:"userSpaceOnUse"};return this.el("pattern",r)};w.use=function(a){return null!=a?(make("use",this.node),a instanceof v&&(a.attr("id")||
a.attr({id:ID()}),a=a.attr("id")),this.el("use",{"xlink:href":a})):v.prototype.use.call(this)};w.text=function(a,d,k){var p={};z(a,"object")?p=a:null!=a&&(p={x:a,y:d,text:k||""});return this.el("text",p)};w.line=function(a,d,k,p){var b={};z(a,"object")?b=a:null!=a&&(b={x1:a,x2:k,y1:d,y2:p});return this.el("line",b)};w.polyline=function(a){1<arguments.length&&(a=Array.prototype.slice.call(arguments,0));var d={};z(a,"object")&&!z(a,"array")?d=a:null!=a&&(d={points:a});return this.el("polyline",d)};
w.polygon=function(a){1<arguments.length&&(a=Array.prototype.slice.call(arguments,0));var d={};z(a,"object")&&!z(a,"array")?d=a:null!=a&&(d={points:a});return this.el("polygon",d)};(function(){function d(){return this.selectAll("stop")}function n(b,d){var f=e("stop"),k={offset:+d+"%"};b=a.color(b);k["stop-color"]=b.hex;1>b.opacity&&(k["stop-opacity"]=b.opacity);e(f,k);this.node.appendChild(f);return this}function u(){if("linearGradient"==this.type){var b=e(this.node,"x1")||0,d=e(this.node,"x2")||
1,f=e(this.node,"y1")||0,k=e(this.node,"y2")||0;return a._.box(b,f,math.abs(d-b),math.abs(k-f))}b=this.node.r||0;return a._.box((this.node.cx||0.5)-b,(this.node.cy||0.5)-b,2*b,2*b)}function p(a,d){function f(a,b){for(var d=(b-u)/(a-w),e=w;e<a;e++)h[e].offset=+(+u+d*(e-w)).toFixed(2);w=a;u=b}var n=k("snap.util.grad.parse",null,d).firstDefined(),p;if(!n)return null;n.params.unshift(a);p="l"==n.type.toLowerCase()?b.apply(0,n.params):q.apply(0,n.params);n.type!=n.type.toLowerCase()&&e(p.node,{gradientUnits:"userSpaceOnUse"});
var h=n.stops,n=h.length,u=0,w=0;n--;for(var v=0;v<n;v++)"offset"in h[v]&&f(v,h[v].offset);h[n].offset=h[n].offset||100;f(n,h[n].offset);for(v=0;v<=n;v++){var y=h[v];p.addStop(y.color,y.offset)}return p}function b(b,k,p,q,w){b=a._.make("linearGradient",b);b.stops=d;b.addStop=n;b.getBBox=u;null!=k&&e(b.node,{x1:k,y1:p,x2:q,y2:w});return b}function q(b,k,p,q,w,h){b=a._.make("radialGradient",b);b.stops=d;b.addStop=n;b.getBBox=u;null!=k&&e(b.node,{cx:k,cy:p,r:q});null!=w&&null!=h&&e(b.node,{fx:w,fy:h});
return b}var e=a._.$;w.gradient=function(a){return p(this.defs,a)};w.gradientLinear=function(a,d,e,f){return b(this.defs,a,d,e,f)};w.gradientRadial=function(a,b,d,e,f){return q(this.defs,a,b,d,e,f)};w.toString=function(){var b=this.node.ownerDocument,d=b.createDocumentFragment(),b=b.createElement("div"),e=this.node.cloneNode(!0);d.appendChild(b);b.appendChild(e);a._.$(e,{xmlns:"http://www.w3.org/2000/svg"});b=b.innerHTML;d.removeChild(d.firstChild);return b};w.clear=function(){for(var a=this.node.firstChild,
b;a;)b=a.nextSibling,"defs"!=a.tagName?a.parentNode.removeChild(a):w.clear.call({node:a}),a=b}})()});C.plugin(function(a,k,y,M){function A(a){var b=A.ps=A.ps||{};b[a]?b[a].sleep=100:b[a]={sleep:100};setTimeout(function(){for(var d in b)b[L](d)&&d!=a&&(b[d].sleep--,!b[d].sleep&&delete b[d])});return b[a]}function w(a,b,d,e){null==a&&(a=b=d=e=0);null==b&&(b=a.y,d=a.width,e=a.height,a=a.x);return{x:a,y:b,width:d,w:d,height:e,h:e,x2:a+d,y2:b+e,cx:a+d/2,cy:b+e/2,r1:F.min(d,e)/2,r2:F.max(d,e)/2,r0:F.sqrt(d*
d+e*e)/2,path:s(a,b,d,e),vb:[a,b,d,e].join(" ")}}function z(){return this.join(",").replace(N,"$1")}function d(a){a=C(a);a.toString=z;return a}function f(a,b,d,h,f,k,l,n,p){if(null==p)return e(a,b,d,h,f,k,l,n);if(0>p||e(a,b,d,h,f,k,l,n)<p)p=void 0;else{var q=0.5,O=1-q,s;for(s=e(a,b,d,h,f,k,l,n,O);0.01<Z(s-p);)q/=2,O+=(s<p?1:-1)*q,s=e(a,b,d,h,f,k,l,n,O);p=O}return u(a,b,d,h,f,k,l,n,p)}function n(b,d){function e(a){return+(+a).toFixed(3)}return a._.cacher(function(a,h,l){a instanceof k&&(a=a.attr("d"));
a=I(a);for(var n,p,D,q,O="",s={},c=0,t=0,r=a.length;t<r;t++){D=a[t];if("M"==D[0])n=+D[1],p=+D[2];else{q=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6]);if(c+q>h){if(d&&!s.start){n=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6],h-c);O+=["C"+e(n.start.x),e(n.start.y),e(n.m.x),e(n.m.y),e(n.x),e(n.y)];if(l)return O;s.start=O;O=["M"+e(n.x),e(n.y)+"C"+e(n.n.x),e(n.n.y),e(n.end.x),e(n.end.y),e(D[5]),e(D[6])].join();c+=q;n=+D[5];p=+D[6];continue}if(!b&&!d)return n=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6],h-c)}c+=q;n=+D[5];p=+D[6]}O+=
D.shift()+D}s.end=O;return n=b?c:d?s:u(n,p,D[0],D[1],D[2],D[3],D[4],D[5],1)},null,a._.clone)}function u(a,b,d,e,h,f,k,l,n){var p=1-n,q=ma(p,3),s=ma(p,2),c=n*n,t=c*n,r=q*a+3*s*n*d+3*p*n*n*h+t*k,q=q*b+3*s*n*e+3*p*n*n*f+t*l,s=a+2*n*(d-a)+c*(h-2*d+a),t=b+2*n*(e-b)+c*(f-2*e+b),x=d+2*n*(h-d)+c*(k-2*h+d),c=e+2*n*(f-e)+c*(l-2*f+e);a=p*a+n*d;b=p*b+n*e;h=p*h+n*k;f=p*f+n*l;l=90-180*F.atan2(s-x,t-c)/S;return{x:r,y:q,m:{x:s,y:t},n:{x:x,y:c},start:{x:a,y:b},end:{x:h,y:f},alpha:l}}function p(b,d,e,h,f,n,k,l){a.is(b,
"array")||(b=[b,d,e,h,f,n,k,l]);b=U.apply(null,b);return w(b.min.x,b.min.y,b.max.x-b.min.x,b.max.y-b.min.y)}function b(a,b,d){return b>=a.x&&b<=a.x+a.width&&d>=a.y&&d<=a.y+a.height}function q(a,d){a=w(a);d=w(d);return b(d,a.x,a.y)||b(d,a.x2,a.y)||b(d,a.x,a.y2)||b(d,a.x2,a.y2)||b(a,d.x,d.y)||b(a,d.x2,d.y)||b(a,d.x,d.y2)||b(a,d.x2,d.y2)||(a.x<d.x2&&a.x>d.x||d.x<a.x2&&d.x>a.x)&&(a.y<d.y2&&a.y>d.y||d.y<a.y2&&d.y>a.y)}function e(a,b,d,e,h,f,n,k,l){null==l&&(l=1);l=(1<l?1:0>l?0:l)/2;for(var p=[-0.1252,
0.1252,-0.3678,0.3678,-0.5873,0.5873,-0.7699,0.7699,-0.9041,0.9041,-0.9816,0.9816],q=[0.2491,0.2491,0.2335,0.2335,0.2032,0.2032,0.1601,0.1601,0.1069,0.1069,0.0472,0.0472],s=0,c=0;12>c;c++)var t=l*p[c]+l,r=t*(t*(-3*a+9*d-9*h+3*n)+6*a-12*d+6*h)-3*a+3*d,t=t*(t*(-3*b+9*e-9*f+3*k)+6*b-12*e+6*f)-3*b+3*e,s=s+q[c]*F.sqrt(r*r+t*t);return l*s}function l(a,b,d){a=I(a);b=I(b);for(var h,f,l,n,k,s,r,O,x,c,t=d?0:[],w=0,v=a.length;w<v;w++)if(x=a[w],"M"==x[0])h=k=x[1],f=s=x[2];else{"C"==x[0]?(x=[h,f].concat(x.slice(1)),
h=x[6],f=x[7]):(x=[h,f,h,f,k,s,k,s],h=k,f=s);for(var G=0,y=b.length;G<y;G++)if(c=b[G],"M"==c[0])l=r=c[1],n=O=c[2];else{"C"==c[0]?(c=[l,n].concat(c.slice(1)),l=c[6],n=c[7]):(c=[l,n,l,n,r,O,r,O],l=r,n=O);var z;var K=x,B=c;z=d;var H=p(K),J=p(B);if(q(H,J)){for(var H=e.apply(0,K),J=e.apply(0,B),H=~~(H/8),J=~~(J/8),U=[],A=[],F={},M=z?0:[],P=0;P<H+1;P++){var C=u.apply(0,K.concat(P/H));U.push({x:C.x,y:C.y,t:P/H})}for(P=0;P<J+1;P++)C=u.apply(0,B.concat(P/J)),A.push({x:C.x,y:C.y,t:P/J});for(P=0;P<H;P++)for(K=
0;K<J;K++){var Q=U[P],L=U[P+1],B=A[K],C=A[K+1],N=0.001>Z(L.x-Q.x)?"y":"x",S=0.001>Z(C.x-B.x)?"y":"x",R;R=Q.x;var Y=Q.y,V=L.x,ea=L.y,fa=B.x,ga=B.y,ha=C.x,ia=C.y;if(W(R,V)<X(fa,ha)||X(R,V)>W(fa,ha)||W(Y,ea)<X(ga,ia)||X(Y,ea)>W(ga,ia))R=void 0;else{var $=(R*ea-Y*V)*(fa-ha)-(R-V)*(fa*ia-ga*ha),aa=(R*ea-Y*V)*(ga-ia)-(Y-ea)*(fa*ia-ga*ha),ja=(R-V)*(ga-ia)-(Y-ea)*(fa-ha);if(ja){var $=$/ja,aa=aa/ja,ja=+$.toFixed(2),ba=+aa.toFixed(2);R=ja<+X(R,V).toFixed(2)||ja>+W(R,V).toFixed(2)||ja<+X(fa,ha).toFixed(2)||
ja>+W(fa,ha).toFixed(2)||ba<+X(Y,ea).toFixed(2)||ba>+W(Y,ea).toFixed(2)||ba<+X(ga,ia).toFixed(2)||ba>+W(ga,ia).toFixed(2)?void 0:{x:$,y:aa}}else R=void 0}R&&F[R.x.toFixed(4)]!=R.y.toFixed(4)&&(F[R.x.toFixed(4)]=R.y.toFixed(4),Q=Q.t+Z((R[N]-Q[N])/(L[N]-Q[N]))*(L.t-Q.t),B=B.t+Z((R[S]-B[S])/(C[S]-B[S]))*(C.t-B.t),0<=Q&&1>=Q&&0<=B&&1>=B&&(z?M++:M.push({x:R.x,y:R.y,t1:Q,t2:B})))}z=M}else z=z?0:[];if(d)t+=z;else{H=0;for(J=z.length;H<J;H++)z[H].segment1=w,z[H].segment2=G,z[H].bez1=x,z[H].bez2=c;t=t.concat(z)}}}return t}
function r(a){var b=A(a);if(b.bbox)return C(b.bbox);if(!a)return w();a=I(a);for(var d=0,e=0,h=[],f=[],l,n=0,k=a.length;n<k;n++)l=a[n],"M"==l[0]?(d=l[1],e=l[2],h.push(d),f.push(e)):(d=U(d,e,l[1],l[2],l[3],l[4],l[5],l[6]),h=h.concat(d.min.x,d.max.x),f=f.concat(d.min.y,d.max.y),d=l[5],e=l[6]);a=X.apply(0,h);l=X.apply(0,f);h=W.apply(0,h);f=W.apply(0,f);f=w(a,l,h-a,f-l);b.bbox=C(f);return f}function s(a,b,d,e,h){if(h)return[["M",+a+ +h,b],["l",d-2*h,0],["a",h,h,0,0,1,h,h],["l",0,e-2*h],["a",h,h,0,0,1,
-h,h],["l",2*h-d,0],["a",h,h,0,0,1,-h,-h],["l",0,2*h-e],["a",h,h,0,0,1,h,-h],["z"] ];a=[["M",a,b],["l",d,0],["l",0,e],["l",-d,0],["z"] ];a.toString=z;return a}function x(a,b,d,e,h){null==h&&null==e&&(e=d);a=+a;b=+b;d=+d;e=+e;if(null!=h){var f=Math.PI/180,l=a+d*Math.cos(-e*f);a+=d*Math.cos(-h*f);var n=b+d*Math.sin(-e*f);b+=d*Math.sin(-h*f);d=[["M",l,n],["A",d,d,0,+(180<h-e),0,a,b] ]}else d=[["M",a,b],["m",0,-e],["a",d,e,0,1,1,0,2*e],["a",d,e,0,1,1,0,-2*e],["z"] ];d.toString=z;return d}function G(b){var e=
A(b);if(e.abs)return d(e.abs);Q(b,"array")&&Q(b&&b[0],"array")||(b=a.parsePathString(b));if(!b||!b.length)return[["M",0,0] ];var h=[],f=0,l=0,n=0,k=0,p=0;"M"==b[0][0]&&(f=+b[0][1],l=+b[0][2],n=f,k=l,p++,h[0]=["M",f,l]);for(var q=3==b.length&&"M"==b[0][0]&&"R"==b[1][0].toUpperCase()&&"Z"==b[2][0].toUpperCase(),s,r,w=p,c=b.length;w<c;w++){h.push(s=[]);r=b[w];p=r[0];if(p!=p.toUpperCase())switch(s[0]=p.toUpperCase(),s[0]){case "A":s[1]=r[1];s[2]=r[2];s[3]=r[3];s[4]=r[4];s[5]=r[5];s[6]=+r[6]+f;s[7]=+r[7]+
l;break;case "V":s[1]=+r[1]+l;break;case "H":s[1]=+r[1]+f;break;case "R":for(var t=[f,l].concat(r.slice(1)),u=2,v=t.length;u<v;u++)t[u]=+t[u]+f,t[++u]=+t[u]+l;h.pop();h=h.concat(P(t,q));break;case "O":h.pop();t=x(f,l,r[1],r[2]);t.push(t[0]);h=h.concat(t);break;case "U":h.pop();h=h.concat(x(f,l,r[1],r[2],r[3]));s=["U"].concat(h[h.length-1].slice(-2));break;case "M":n=+r[1]+f,k=+r[2]+l;default:for(u=1,v=r.length;u<v;u++)s[u]=+r[u]+(u%2?f:l)}else if("R"==p)t=[f,l].concat(r.slice(1)),h.pop(),h=h.concat(P(t,
q)),s=["R"].concat(r.slice(-2));else if("O"==p)h.pop(),t=x(f,l,r[1],r[2]),t.push(t[0]),h=h.concat(t);else if("U"==p)h.pop(),h=h.concat(x(f,l,r[1],r[2],r[3])),s=["U"].concat(h[h.length-1].slice(-2));else for(t=0,u=r.length;t<u;t++)s[t]=r[t];p=p.toUpperCase();if("O"!=p)switch(s[0]){case "Z":f=+n;l=+k;break;case "H":f=s[1];break;case "V":l=s[1];break;case "M":n=s[s.length-2],k=s[s.length-1];default:f=s[s.length-2],l=s[s.length-1]}}h.toString=z;e.abs=d(h);return h}function h(a,b,d,e){return[a,b,d,e,d,
e]}function J(a,b,d,e,h,f){var l=1/3,n=2/3;return[l*a+n*d,l*b+n*e,l*h+n*d,l*f+n*e,h,f]}function K(b,d,e,h,f,l,n,k,p,s){var r=120*S/180,q=S/180*(+f||0),c=[],t,x=a._.cacher(function(a,b,c){var d=a*F.cos(c)-b*F.sin(c);a=a*F.sin(c)+b*F.cos(c);return{x:d,y:a}});if(s)v=s[0],t=s[1],l=s[2],u=s[3];else{t=x(b,d,-q);b=t.x;d=t.y;t=x(k,p,-q);k=t.x;p=t.y;F.cos(S/180*f);F.sin(S/180*f);t=(b-k)/2;v=(d-p)/2;u=t*t/(e*e)+v*v/(h*h);1<u&&(u=F.sqrt(u),e*=u,h*=u);var u=e*e,w=h*h,u=(l==n?-1:1)*F.sqrt(Z((u*w-u*v*v-w*t*t)/
(u*v*v+w*t*t)));l=u*e*v/h+(b+k)/2;var u=u*-h*t/e+(d+p)/2,v=F.asin(((d-u)/h).toFixed(9));t=F.asin(((p-u)/h).toFixed(9));v=b<l?S-v:v;t=k<l?S-t:t;0>v&&(v=2*S+v);0>t&&(t=2*S+t);n&&v>t&&(v-=2*S);!n&&t>v&&(t-=2*S)}if(Z(t-v)>r){var c=t,w=k,G=p;t=v+r*(n&&t>v?1:-1);k=l+e*F.cos(t);p=u+h*F.sin(t);c=K(k,p,e,h,f,0,n,w,G,[t,c,l,u])}l=t-v;f=F.cos(v);r=F.sin(v);n=F.cos(t);t=F.sin(t);l=F.tan(l/4);e=4/3*e*l;l*=4/3*h;h=[b,d];b=[b+e*r,d-l*f];d=[k+e*t,p-l*n];k=[k,p];b[0]=2*h[0]-b[0];b[1]=2*h[1]-b[1];if(s)return[b,d,k].concat(c);
c=[b,d,k].concat(c).join().split(",");s=[];k=0;for(p=c.length;k<p;k++)s[k]=k%2?x(c[k-1],c[k],q).y:x(c[k],c[k+1],q).x;return s}function U(a,b,d,e,h,f,l,k){for(var n=[],p=[[],[] ],s,r,c,t,q=0;2>q;++q)0==q?(r=6*a-12*d+6*h,s=-3*a+9*d-9*h+3*l,c=3*d-3*a):(r=6*b-12*e+6*f,s=-3*b+9*e-9*f+3*k,c=3*e-3*b),1E-12>Z(s)?1E-12>Z(r)||(s=-c/r,0<s&&1>s&&n.push(s)):(t=r*r-4*c*s,c=F.sqrt(t),0>t||(t=(-r+c)/(2*s),0<t&&1>t&&n.push(t),s=(-r-c)/(2*s),0<s&&1>s&&n.push(s)));for(r=q=n.length;q--;)s=n[q],c=1-s,p[0][q]=c*c*c*a+3*
c*c*s*d+3*c*s*s*h+s*s*s*l,p[1][q]=c*c*c*b+3*c*c*s*e+3*c*s*s*f+s*s*s*k;p[0][r]=a;p[1][r]=b;p[0][r+1]=l;p[1][r+1]=k;p[0].length=p[1].length=r+2;return{min:{x:X.apply(0,p[0]),y:X.apply(0,p[1])},max:{x:W.apply(0,p[0]),y:W.apply(0,p[1])}}}function I(a,b){var e=!b&&A(a);if(!b&&e.curve)return d(e.curve);var f=G(a),l=b&&G(b),n={x:0,y:0,bx:0,by:0,X:0,Y:0,qx:null,qy:null},k={x:0,y:0,bx:0,by:0,X:0,Y:0,qx:null,qy:null},p=function(a,b,c){if(!a)return["C",b.x,b.y,b.x,b.y,b.x,b.y];a[0]in{T:1,Q:1}||(b.qx=b.qy=null);
switch(a[0]){case "M":b.X=a[1];b.Y=a[2];break;case "A":a=["C"].concat(K.apply(0,[b.x,b.y].concat(a.slice(1))));break;case "S":"C"==c||"S"==c?(c=2*b.x-b.bx,b=2*b.y-b.by):(c=b.x,b=b.y);a=["C",c,b].concat(a.slice(1));break;case "T":"Q"==c||"T"==c?(b.qx=2*b.x-b.qx,b.qy=2*b.y-b.qy):(b.qx=b.x,b.qy=b.y);a=["C"].concat(J(b.x,b.y,b.qx,b.qy,a[1],a[2]));break;case "Q":b.qx=a[1];b.qy=a[2];a=["C"].concat(J(b.x,b.y,a[1],a[2],a[3],a[4]));break;case "L":a=["C"].concat(h(b.x,b.y,a[1],a[2]));break;case "H":a=["C"].concat(h(b.x,
b.y,a[1],b.y));break;case "V":a=["C"].concat(h(b.x,b.y,b.x,a[1]));break;case "Z":a=["C"].concat(h(b.x,b.y,b.X,b.Y))}return a},s=function(a,b){if(7<a[b].length){a[b].shift();for(var c=a[b];c.length;)q[b]="A",l&&(u[b]="A"),a.splice(b++,0,["C"].concat(c.splice(0,6)));a.splice(b,1);v=W(f.length,l&&l.length||0)}},r=function(a,b,c,d,e){a&&b&&"M"==a[e][0]&&"M"!=b[e][0]&&(b.splice(e,0,["M",d.x,d.y]),c.bx=0,c.by=0,c.x=a[e][1],c.y=a[e][2],v=W(f.length,l&&l.length||0))},q=[],u=[],c="",t="",x=0,v=W(f.length,
l&&l.length||0);for(;x<v;x++){f[x]&&(c=f[x][0]);"C"!=c&&(q[x]=c,x&&(t=q[x-1]));f[x]=p(f[x],n,t);"A"!=q[x]&&"C"==c&&(q[x]="C");s(f,x);l&&(l[x]&&(c=l[x][0]),"C"!=c&&(u[x]=c,x&&(t=u[x-1])),l[x]=p(l[x],k,t),"A"!=u[x]&&"C"==c&&(u[x]="C"),s(l,x));r(f,l,n,k,x);r(l,f,k,n,x);var w=f[x],z=l&&l[x],y=w.length,U=l&&z.length;n.x=w[y-2];n.y=w[y-1];n.bx=$(w[y-4])||n.x;n.by=$(w[y-3])||n.y;k.bx=l&&($(z[U-4])||k.x);k.by=l&&($(z[U-3])||k.y);k.x=l&&z[U-2];k.y=l&&z[U-1]}l||(e.curve=d(f));return l?[f,l]:f}function P(a,
b){for(var d=[],e=0,h=a.length;h-2*!b>e;e+=2){var f=[{x:+a[e-2],y:+a[e-1]},{x:+a[e],y:+a[e+1]},{x:+a[e+2],y:+a[e+3]},{x:+a[e+4],y:+a[e+5]}];b?e?h-4==e?f[3]={x:+a[0],y:+a[1]}:h-2==e&&(f[2]={x:+a[0],y:+a[1]},f[3]={x:+a[2],y:+a[3]}):f[0]={x:+a[h-2],y:+a[h-1]}:h-4==e?f[3]=f[2]:e||(f[0]={x:+a[e],y:+a[e+1]});d.push(["C",(-f[0].x+6*f[1].x+f[2].x)/6,(-f[0].y+6*f[1].y+f[2].y)/6,(f[1].x+6*f[2].x-f[3].x)/6,(f[1].y+6*f[2].y-f[3].y)/6,f[2].x,f[2].y])}return d}y=k.prototype;var Q=a.is,C=a._.clone,L="hasOwnProperty",
N=/,?([a-z]),?/gi,$=parseFloat,F=Math,S=F.PI,X=F.min,W=F.max,ma=F.pow,Z=F.abs;M=n(1);var na=n(),ba=n(0,1),V=a._unit2px;a.path=A;a.path.getTotalLength=M;a.path.getPointAtLength=na;a.path.getSubpath=function(a,b,d){if(1E-6>this.getTotalLength(a)-d)return ba(a,b).end;a=ba(a,d,1);return b?ba(a,b).end:a};y.getTotalLength=function(){if(this.node.getTotalLength)return this.node.getTotalLength()};y.getPointAtLength=function(a){return na(this.attr("d"),a)};y.getSubpath=function(b,d){return a.path.getSubpath(this.attr("d"),
b,d)};a._.box=w;a.path.findDotsAtSegment=u;a.path.bezierBBox=p;a.path.isPointInsideBBox=b;a.path.isBBoxIntersect=q;a.path.intersection=function(a,b){return l(a,b)};a.path.intersectionNumber=function(a,b){return l(a,b,1)};a.path.isPointInside=function(a,d,e){var h=r(a);return b(h,d,e)&&1==l(a,[["M",d,e],["H",h.x2+10] ],1)%2};a.path.getBBox=r;a.path.get={path:function(a){return a.attr("path")},circle:function(a){a=V(a);return x(a.cx,a.cy,a.r)},ellipse:function(a){a=V(a);return x(a.cx||0,a.cy||0,a.rx,
a.ry)},rect:function(a){a=V(a);return s(a.x||0,a.y||0,a.width,a.height,a.rx,a.ry)},image:function(a){a=V(a);return s(a.x||0,a.y||0,a.width,a.height)},line:function(a){return"M"+[a.attr("x1")||0,a.attr("y1")||0,a.attr("x2"),a.attr("y2")]},polyline:function(a){return"M"+a.attr("points")},polygon:function(a){return"M"+a.attr("points")+"z"},deflt:function(a){a=a.node.getBBox();return s(a.x,a.y,a.width,a.height)}};a.path.toRelative=function(b){var e=A(b),h=String.prototype.toLowerCase;if(e.rel)return d(e.rel);
a.is(b,"array")&&a.is(b&&b[0],"array")||(b=a.parsePathString(b));var f=[],l=0,n=0,k=0,p=0,s=0;"M"==b[0][0]&&(l=b[0][1],n=b[0][2],k=l,p=n,s++,f.push(["M",l,n]));for(var r=b.length;s<r;s++){var q=f[s]=[],x=b[s];if(x[0]!=h.call(x[0]))switch(q[0]=h.call(x[0]),q[0]){case "a":q[1]=x[1];q[2]=x[2];q[3]=x[3];q[4]=x[4];q[5]=x[5];q[6]=+(x[6]-l).toFixed(3);q[7]=+(x[7]-n).toFixed(3);break;case "v":q[1]=+(x[1]-n).toFixed(3);break;case "m":k=x[1],p=x[2];default:for(var c=1,t=x.length;c<t;c++)q[c]=+(x[c]-(c%2?l:
n)).toFixed(3)}else for(f[s]=[],"m"==x[0]&&(k=x[1]+l,p=x[2]+n),q=0,c=x.length;q<c;q++)f[s][q]=x[q];x=f[s].length;switch(f[s][0]){case "z":l=k;n=p;break;case "h":l+=+f[s][x-1];break;case "v":n+=+f[s][x-1];break;default:l+=+f[s][x-2],n+=+f[s][x-1]}}f.toString=z;e.rel=d(f);return f};a.path.toAbsolute=G;a.path.toCubic=I;a.path.map=function(a,b){if(!b)return a;var d,e,h,f,l,n,k;a=I(a);h=0;for(l=a.length;h<l;h++)for(k=a[h],f=1,n=k.length;f<n;f+=2)d=b.x(k[f],k[f+1]),e=b.y(k[f],k[f+1]),k[f]=d,k[f+1]=e;return a};
a.path.toString=z;a.path.clone=d});C.plugin(function(a,v,y,C){var A=Math.max,w=Math.min,z=function(a){this.items=[];this.bindings={};this.length=0;this.type="set";if(a)for(var f=0,n=a.length;f<n;f++)a[f]&&(this[this.items.length]=this.items[this.items.length]=a[f],this.length++)};v=z.prototype;v.push=function(){for(var a,f,n=0,k=arguments.length;n<k;n++)if(a=arguments[n])f=this.items.length,this[f]=this.items[f]=a,this.length++;return this};v.pop=function(){this.length&&delete this[this.length--];
return this.items.pop()};v.forEach=function(a,f){for(var n=0,k=this.items.length;n<k&&!1!==a.call(f,this.items[n],n);n++);return this};v.animate=function(d,f,n,u){"function"!=typeof n||n.length||(u=n,n=L.linear);d instanceof a._.Animation&&(u=d.callback,n=d.easing,f=n.dur,d=d.attr);var p=arguments;if(a.is(d,"array")&&a.is(p[p.length-1],"array"))var b=!0;var q,e=function(){q?this.b=q:q=this.b},l=0,r=u&&function(){l++==this.length&&u.call(this)};return this.forEach(function(a,l){k.once("snap.animcreated."+
a.id,e);b?p[l]&&a.animate.apply(a,p[l]):a.animate(d,f,n,r)})};v.remove=function(){for(;this.length;)this.pop().remove();return this};v.bind=function(a,f,k){var u={};if("function"==typeof f)this.bindings[a]=f;else{var p=k||a;this.bindings[a]=function(a){u[p]=a;f.attr(u)}}return this};v.attr=function(a){var f={},k;for(k in a)if(this.bindings[k])this.bindings[k](a[k]);else f[k]=a[k];a=0;for(k=this.items.length;a<k;a++)this.items[a].attr(f);return this};v.clear=function(){for(;this.length;)this.pop()};
v.splice=function(a,f,k){a=0>a?A(this.length+a,0):a;f=A(0,w(this.length-a,f));var u=[],p=[],b=[],q;for(q=2;q<arguments.length;q++)b.push(arguments[q]);for(q=0;q<f;q++)p.push(this[a+q]);for(;q<this.length-a;q++)u.push(this[a+q]);var e=b.length;for(q=0;q<e+u.length;q++)this.items[a+q]=this[a+q]=q<e?b[q]:u[q-e];for(q=this.items.length=this.length-=f-e;this[q];)delete this[q++];return new z(p)};v.exclude=function(a){for(var f=0,k=this.length;f<k;f++)if(this[f]==a)return this.splice(f,1),!0;return!1};
v.insertAfter=function(a){for(var f=this.items.length;f--;)this.items[f].insertAfter(a);return this};v.getBBox=function(){for(var a=[],f=[],k=[],u=[],p=this.items.length;p--;)if(!this.items[p].removed){var b=this.items[p].getBBox();a.push(b.x);f.push(b.y);k.push(b.x+b.width);u.push(b.y+b.height)}a=w.apply(0,a);f=w.apply(0,f);k=A.apply(0,k);u=A.apply(0,u);return{x:a,y:f,x2:k,y2:u,width:k-a,height:u-f,cx:a+(k-a)/2,cy:f+(u-f)/2}};v.clone=function(a){a=new z;for(var f=0,k=this.items.length;f<k;f++)a.push(this.items[f].clone());
return a};v.toString=function(){return"Snap\u2018s set"};v.type="set";a.set=function(){var a=new z;arguments.length&&a.push.apply(a,Array.prototype.slice.call(arguments,0));return a}});C.plugin(function(a,v,y,C){function A(a){var b=a[0];switch(b.toLowerCase()){case "t":return[b,0,0];case "m":return[b,1,0,0,1,0,0];case "r":return 4==a.length?[b,0,a[2],a[3] ]:[b,0];case "s":return 5==a.length?[b,1,1,a[3],a[4] ]:3==a.length?[b,1,1]:[b,1]}}function w(b,d,f){d=q(d).replace(/\.{3}|\u2026/g,b);b=a.parseTransformString(b)||
[];d=a.parseTransformString(d)||[];for(var k=Math.max(b.length,d.length),p=[],v=[],h=0,w,z,y,I;h<k;h++){y=b[h]||A(d[h]);I=d[h]||A(y);if(y[0]!=I[0]||"r"==y[0].toLowerCase()&&(y[2]!=I[2]||y[3]!=I[3])||"s"==y[0].toLowerCase()&&(y[3]!=I[3]||y[4]!=I[4])){b=a._.transform2matrix(b,f());d=a._.transform2matrix(d,f());p=[["m",b.a,b.b,b.c,b.d,b.e,b.f] ];v=[["m",d.a,d.b,d.c,d.d,d.e,d.f] ];break}p[h]=[];v[h]=[];w=0;for(z=Math.max(y.length,I.length);w<z;w++)w in y&&(p[h][w]=y[w]),w in I&&(v[h][w]=I[w])}return{from:u(p),
to:u(v),f:n(p)}}function z(a){return a}function d(a){return function(b){return+b.toFixed(3)+a}}function f(b){return a.rgb(b[0],b[1],b[2])}function n(a){var b=0,d,f,k,n,h,p,q=[];d=0;for(f=a.length;d<f;d++){h="[";p=['"'+a[d][0]+'"'];k=1;for(n=a[d].length;k<n;k++)p[k]="val["+b++ +"]";h+=p+"]";q[d]=h}return Function("val","return Snap.path.toString.call(["+q+"])")}function u(a){for(var b=[],d=0,f=a.length;d<f;d++)for(var k=1,n=a[d].length;k<n;k++)b.push(a[d][k]);return b}var p={},b=/[a-z]+$/i,q=String;
p.stroke=p.fill="colour";v.prototype.equal=function(a,b){return k("snap.util.equal",this,a,b).firstDefined()};k.on("snap.util.equal",function(e,k){var r,s;r=q(this.attr(e)||"");var x=this;if(r==+r&&k==+k)return{from:+r,to:+k,f:z};if("colour"==p[e])return r=a.color(r),s=a.color(k),{from:[r.r,r.g,r.b,r.opacity],to:[s.r,s.g,s.b,s.opacity],f:f};if("transform"==e||"gradientTransform"==e||"patternTransform"==e)return k instanceof a.Matrix&&(k=k.toTransformString()),a._.rgTransform.test(k)||(k=a._.svgTransform2string(k)),
w(r,k,function(){return x.getBBox(1)});if("d"==e||"path"==e)return r=a.path.toCubic(r,k),{from:u(r[0]),to:u(r[1]),f:n(r[0])};if("points"==e)return r=q(r).split(a._.separator),s=q(k).split(a._.separator),{from:r,to:s,f:function(a){return a}};aUnit=r.match(b);s=q(k).match(b);return aUnit&&aUnit==s?{from:parseFloat(r),to:parseFloat(k),f:d(aUnit)}:{from:this.asPX(e),to:this.asPX(e,k),f:z}})});C.plugin(function(a,v,y,C){var A=v.prototype,w="createTouch"in C.doc;v="click dblclick mousedown mousemove mouseout mouseover mouseup touchstart touchmove touchend touchcancel".split(" ");
var z={mousedown:"touchstart",mousemove:"touchmove",mouseup:"touchend"},d=function(a,b){var d="y"==a?"scrollTop":"scrollLeft",e=b&&b.node?b.node.ownerDocument:C.doc;return e[d in e.documentElement?"documentElement":"body"][d]},f=function(){this.returnValue=!1},n=function(){return this.originalEvent.preventDefault()},u=function(){this.cancelBubble=!0},p=function(){return this.originalEvent.stopPropagation()},b=function(){if(C.doc.addEventListener)return function(a,b,e,f){var k=w&&z[b]?z[b]:b,l=function(k){var l=
d("y",f),q=d("x",f);if(w&&z.hasOwnProperty(b))for(var r=0,u=k.targetTouches&&k.targetTouches.length;r<u;r++)if(k.targetTouches[r].target==a||a.contains(k.targetTouches[r].target)){u=k;k=k.targetTouches[r];k.originalEvent=u;k.preventDefault=n;k.stopPropagation=p;break}return e.call(f,k,k.clientX+q,k.clientY+l)};b!==k&&a.addEventListener(b,l,!1);a.addEventListener(k,l,!1);return function(){b!==k&&a.removeEventListener(b,l,!1);a.removeEventListener(k,l,!1);return!0}};if(C.doc.attachEvent)return function(a,
b,e,h){var k=function(a){a=a||h.node.ownerDocument.window.event;var b=d("y",h),k=d("x",h),k=a.clientX+k,b=a.clientY+b;a.preventDefault=a.preventDefault||f;a.stopPropagation=a.stopPropagation||u;return e.call(h,a,k,b)};a.attachEvent("on"+b,k);return function(){a.detachEvent("on"+b,k);return!0}}}(),q=[],e=function(a){for(var b=a.clientX,e=a.clientY,f=d("y"),l=d("x"),n,p=q.length;p--;){n=q[p];if(w)for(var r=a.touches&&a.touches.length,u;r--;){if(u=a.touches[r],u.identifier==n.el._drag.id||n.el.node.contains(u.target)){b=
u.clientX;e=u.clientY;(a.originalEvent?a.originalEvent:a).preventDefault();break}}else a.preventDefault();b+=l;e+=f;k("snap.drag.move."+n.el.id,n.move_scope||n.el,b-n.el._drag.x,e-n.el._drag.y,b,e,a)}},l=function(b){a.unmousemove(e).unmouseup(l);for(var d=q.length,f;d--;)f=q[d],f.el._drag={},k("snap.drag.end."+f.el.id,f.end_scope||f.start_scope||f.move_scope||f.el,b);q=[]};for(y=v.length;y--;)(function(d){a[d]=A[d]=function(e,f){a.is(e,"function")&&(this.events=this.events||[],this.events.push({name:d,
f:e,unbind:b(this.node||document,d,e,f||this)}));return this};a["un"+d]=A["un"+d]=function(a){for(var b=this.events||[],e=b.length;e--;)if(b[e].name==d&&(b[e].f==a||!a)){b[e].unbind();b.splice(e,1);!b.length&&delete this.events;break}return this}})(v[y]);A.hover=function(a,b,d,e){return this.mouseover(a,d).mouseout(b,e||d)};A.unhover=function(a,b){return this.unmouseover(a).unmouseout(b)};var r=[];A.drag=function(b,d,f,h,n,p){function u(r,v,w){(r.originalEvent||r).preventDefault();this._drag.x=v;
this._drag.y=w;this._drag.id=r.identifier;!q.length&&a.mousemove(e).mouseup(l);q.push({el:this,move_scope:h,start_scope:n,end_scope:p});d&&k.on("snap.drag.start."+this.id,d);b&&k.on("snap.drag.move."+this.id,b);f&&k.on("snap.drag.end."+this.id,f);k("snap.drag.start."+this.id,n||h||this,v,w,r)}if(!arguments.length){var v;return this.drag(function(a,b){this.attr({transform:v+(v?"T":"t")+[a,b]})},function(){v=this.transform().local})}this._drag={};r.push({el:this,start:u});this.mousedown(u);return this};
A.undrag=function(){for(var b=r.length;b--;)r[b].el==this&&(this.unmousedown(r[b].start),r.splice(b,1),k.unbind("snap.drag.*."+this.id));!r.length&&a.unmousemove(e).unmouseup(l);return this}});C.plugin(function(a,v,y,C){y=y.prototype;var A=/^\s*url\((.+)\)/,w=String,z=a._.$;a.filter={};y.filter=function(d){var f=this;"svg"!=f.type&&(f=f.paper);d=a.parse(w(d));var k=a._.id(),u=z("filter");z(u,{id:k,filterUnits:"userSpaceOnUse"});u.appendChild(d.node);f.defs.appendChild(u);return new v(u)};k.on("snap.util.getattr.filter",
function(){k.stop();var d=z(this.node,"filter");if(d)return(d=w(d).match(A))&&a.select(d[1])});k.on("snap.util.attr.filter",function(d){if(d instanceof v&&"filter"==d.type){k.stop();var f=d.node.id;f||(z(d.node,{id:d.id}),f=d.id);z(this.node,{filter:a.url(f)})}d&&"none"!=d||(k.stop(),this.node.removeAttribute("filter"))});a.filter.blur=function(d,f){null==d&&(d=2);return a.format('<feGaussianBlur stdDeviation="{def}"/>',{def:null==f?d:[d,f]})};a.filter.blur.toString=function(){return this()};a.filter.shadow=
function(d,f,k,u,p){"string"==typeof k&&(p=u=k,k=4);"string"!=typeof u&&(p=u,u="#000");null==k&&(k=4);null==p&&(p=1);null==d&&(d=0,f=2);null==f&&(f=d);u=a.color(u||"#000");return a.format('<feGaussianBlur in="SourceAlpha" stdDeviation="{blur}"/><feOffset dx="{dx}" dy="{dy}" result="offsetblur"/><feFlood flood-color="{color}"/><feComposite in2="offsetblur" operator="in"/><feComponentTransfer><feFuncA type="linear" slope="{opacity}"/></feComponentTransfer><feMerge><feMergeNode/><feMergeNode in="SourceGraphic"/></feMerge>',
{color:u,dx:d,dy:f,blur:k,opacity:p})};a.filter.shadow.toString=function(){return this()};a.filter.grayscale=function(d){null==d&&(d=1);return a.format('<feColorMatrix type="matrix" values="{a} {b} {c} 0 0 {d} {e} {f} 0 0 {g} {b} {h} 0 0 0 0 0 1 0"/>',{a:0.2126+0.7874*(1-d),b:0.7152-0.7152*(1-d),c:0.0722-0.0722*(1-d),d:0.2126-0.2126*(1-d),e:0.7152+0.2848*(1-d),f:0.0722-0.0722*(1-d),g:0.2126-0.2126*(1-d),h:0.0722+0.9278*(1-d)})};a.filter.grayscale.toString=function(){return this()};a.filter.sepia=
function(d){null==d&&(d=1);return a.format('<feColorMatrix type="matrix" values="{a} {b} {c} 0 0 {d} {e} {f} 0 0 {g} {h} {i} 0 0 0 0 0 1 0"/>',{a:0.393+0.607*(1-d),b:0.769-0.769*(1-d),c:0.189-0.189*(1-d),d:0.349-0.349*(1-d),e:0.686+0.314*(1-d),f:0.168-0.168*(1-d),g:0.272-0.272*(1-d),h:0.534-0.534*(1-d),i:0.131+0.869*(1-d)})};a.filter.sepia.toString=function(){return this()};a.filter.saturate=function(d){null==d&&(d=1);return a.format('<feColorMatrix type="saturate" values="{amount}"/>',{amount:1-
d})};a.filter.saturate.toString=function(){return this()};a.filter.hueRotate=function(d){return a.format('<feColorMatrix type="hueRotate" values="{angle}"/>',{angle:d||0})};a.filter.hueRotate.toString=function(){return this()};a.filter.invert=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="table" tableValues="{amount} {amount2}"/><feFuncG type="table" tableValues="{amount} {amount2}"/><feFuncB type="table" tableValues="{amount} {amount2}"/></feComponentTransfer>',{amount:d,
amount2:1-d})};a.filter.invert.toString=function(){return this()};a.filter.brightness=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="linear" slope="{amount}"/><feFuncG type="linear" slope="{amount}"/><feFuncB type="linear" slope="{amount}"/></feComponentTransfer>',{amount:d})};a.filter.brightness.toString=function(){return this()};a.filter.contrast=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="linear" slope="{amount}" intercept="{amount2}"/><feFuncG type="linear" slope="{amount}" intercept="{amount2}"/><feFuncB type="linear" slope="{amount}" intercept="{amount2}"/></feComponentTransfer>',
{amount:d,amount2:0.5-d/2})};a.filter.contrast.toString=function(){return this()}});return C});

]]> </script>
</svg>


<pre class='language-julia'><code class='language-julia'>mat(qft3_first)</code></pre>
<pre id='var-hash601846' class='code-output documenter-example-output'>8×8 SparseMatrixCSC{ComplexF64, Int64} with 16 stored entries:
 0.707107+0.0im   0.707107+0.0im  …           ⋅               ⋅    
 0.707107+0.0im  -0.707107+0.0im              ⋅               ⋅    
          ⋅                ⋅                  ⋅               ⋅    
          ⋅                ⋅                  ⋅               ⋅    
          ⋅                ⋅                  ⋅               ⋅    
          ⋅                ⋅      …           ⋅               ⋅    
          ⋅                ⋅         0.707107+0.0im  0.707107+0.0im
          ⋅                ⋅             -0.5+0.5im       0.5-0.5im</pre>

<pre class='language-julia'><code class='language-julia'>mat(qft3)</code></pre>
<pre id='var-hash164697' class='code-output documenter-example-output'>8×8 SparseMatrixCSC{ComplexF64, Int64} with 64 stored entries:
 0.353553+0.0im   0.353553+0.0im      0.353553+0.0im       …      0.353553+0.0im
 0.353553+0.0im  -0.353553+0.0im   2.16489e-17+0.353553im             0.25-0.25im
 0.353553+0.0im   0.353553+0.0im     -0.353553+0.0im          -2.16489e-17-0.353553im
 0.353553+0.0im  -0.353553+0.0im  -2.16489e-17-0.353553im            -0.25-0.25im
 0.353553+0.0im   0.353553+0.0im      0.353553+0.0im             -0.353553+0.0im
 0.353553+0.0im  -0.353553+0.0im   2.16489e-17+0.353553im  …         -0.25+0.25im
 0.353553+0.0im   0.353553+0.0im     -0.353553+0.0im           2.16489e-17+0.353553im
 0.353553+0.0im  -0.353553+0.0im  -2.16489e-17-0.353553im             0.25+0.25im</pre>


<div class="markdown"><h3>Use matrix as a block</h3>
</div>

<pre class='language-julia'><code class='language-julia'>random_gate = matblock(rand_unitary(2); tag="random gate")</code></pre>
<pre id='var-random_gate' class='code-output documenter-example-output'>random gate</pre>

<pre class='language-julia'><code class='language-julia'>mat(put(3, 2=&gt;random_gate))</code></pre>
<pre id='var-hash183324' class='code-output documenter-example-output'>8×8 SparseMatrixCSC{ComplexF64, Int64} with 16 stored entries:
 -0.724423-0.483425im            ⋅           …             ⋅    
           ⋅           -0.724423-0.483425im                ⋅    
 0.0647939-0.487148im            ⋅                         ⋅    
           ⋅           0.0647939-0.487148im                ⋅    
           ⋅                     ⋅                         ⋅    
           ⋅                     ⋅           …  -0.0692378-0.486536im
           ⋅                     ⋅                         ⋅    
           ⋅                     ⋅                -0.71998+0.490017im</pre>

<pre class='language-julia'><code class='language-julia'>vizcircuit(put(3, 2=&gt;random_gate))</code></pre>
<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:xlink="http://www.w3.org/1999/xlink"
     version="1.2"
     width="63.75mm" height="56.25mm" viewBox="0 0 63.75 56.25"
     stroke="none"
     fill="#000000"
     stroke-width="0.3"
     font-size="3.88"

     id="img-cfbfb39d">
<defs>
  <marker id="arrow" markerWidth="15" markerHeight="7" refX="5" refY="3.5" orient="auto" markerUnits="strokeWidth">
    <path d="M0,0 L15,3.5 L0,7 z" stroke="context-stroke" fill="context-stroke"/>
  </marker>
</defs>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-cfbfb39d-1">
  <g transform="translate(33.75,25)">
    <path d="M-12.72,-5 L12.72,-5 12.72,5 -12.72,5  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-cfbfb39d-2">
  <g transform="translate(33.75,43.75)">
    <path fill="none" d="M21.25,0 L-21.25,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-cfbfb39d-3">
  <g transform="translate(55,25)">
    <path fill="none" d="M0,0 L0,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-cfbfb39d-4">
  <g transform="translate(33.75,6.25)">
    <path fill="none" d="M21.25,0 L-21.25,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-cfbfb39d-5">
  <g transform="translate(50.73,25)">
    <path fill="none" d="M-4.27,0 L4.27,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-cfbfb39d-6">
  <g transform="translate(16.77,25)">
    <path fill="none" d="M-4.27,0 L4.27,0 " class="primitive"/>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="3.53" id="img-cfbfb39d-7">
  <g transform="translate(33.75,25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">random gate</text>
    </g>
  </g>
</g>
<script> <![CDATA[
(function(N){var k=/[\.\/]/,L=/\s*,\s*/,C=function(a,d){return a-d},a,v,y={n:{}},M=function(){for(var a=0,d=this.length;a<d;a++)if("undefined"!=typeof this[a])return this[a]},A=function(){for(var a=this.length;--a;)if("undefined"!=typeof this[a])return this[a]},w=function(k,d){k=String(k);var f=v,n=Array.prototype.slice.call(arguments,2),u=w.listeners(k),p=0,b,q=[],e={},l=[],r=a;l.firstDefined=M;l.lastDefined=A;a=k;for(var s=v=0,x=u.length;s<x;s++)"zIndex"in u[s]&&(q.push(u[s].zIndex),0>u[s].zIndex&&
(e[u[s].zIndex]=u[s]));for(q.sort(C);0>q[p];)if(b=e[q[p++] ],l.push(b.apply(d,n)),v)return v=f,l;for(s=0;s<x;s++)if(b=u[s],"zIndex"in b)if(b.zIndex==q[p]){l.push(b.apply(d,n));if(v)break;do if(p++,(b=e[q[p] ])&&l.push(b.apply(d,n)),v)break;while(b)}else e[b.zIndex]=b;else if(l.push(b.apply(d,n)),v)break;v=f;a=r;return l};w._events=y;w.listeners=function(a){a=a.split(k);var d=y,f,n,u,p,b,q,e,l=[d],r=[];u=0;for(p=a.length;u<p;u++){e=[];b=0;for(q=l.length;b<q;b++)for(d=l[b].n,f=[d[a[u] ],d["*"] ],n=2;n--;)if(d=
f[n])e.push(d),r=r.concat(d.f||[]);l=e}return r};w.on=function(a,d){a=String(a);if("function"!=typeof d)return function(){};for(var f=a.split(L),n=0,u=f.length;n<u;n++)(function(a){a=a.split(k);for(var b=y,f,e=0,l=a.length;e<l;e++)b=b.n,b=b.hasOwnProperty(a[e])&&b[a[e] ]||(b[a[e] ]={n:{}});b.f=b.f||[];e=0;for(l=b.f.length;e<l;e++)if(b.f[e]==d){f=!0;break}!f&&b.f.push(d)})(f[n]);return function(a){+a==+a&&(d.zIndex=+a)}};w.f=function(a){var d=[].slice.call(arguments,1);return function(){w.apply(null,
[a,null].concat(d).concat([].slice.call(arguments,0)))}};w.stop=function(){v=1};w.nt=function(k){return k?(new RegExp("(?:\\.|\\/|^)"+k+"(?:\\.|\\/|$)")).test(a):a};w.nts=function(){return a.split(k)};w.off=w.unbind=function(a,d){if(a){var f=a.split(L);if(1<f.length)for(var n=0,u=f.length;n<u;n++)w.off(f[n],d);else{for(var f=a.split(k),p,b,q,e,l=[y],n=0,u=f.length;n<u;n++)for(e=0;e<l.length;e+=q.length-2){q=[e,1];p=l[e].n;if("*"!=f[n])p[f[n] ]&&q.push(p[f[n] ]);else for(b in p)p.hasOwnProperty(b)&&
q.push(p[b]);l.splice.apply(l,q)}n=0;for(u=l.length;n<u;n++)for(p=l[n];p.n;){if(d){if(p.f){e=0;for(f=p.f.length;e<f;e++)if(p.f[e]==d){p.f.splice(e,1);break}!p.f.length&&delete p.f}for(b in p.n)if(p.n.hasOwnProperty(b)&&p.n[b].f){q=p.n[b].f;e=0;for(f=q.length;e<f;e++)if(q[e]==d){q.splice(e,1);break}!q.length&&delete p.n[b].f}}else for(b in delete p.f,p.n)p.n.hasOwnProperty(b)&&p.n[b].f&&delete p.n[b].f;p=p.n}}}else w._events=y={n:{}}};w.once=function(a,d){var f=function(){w.unbind(a,f);return d.apply(this,
arguments)};return w.on(a,f)};w.version="0.4.2";w.toString=function(){return"You are running Eve 0.4.2"};"undefined"!=typeof module&&module.exports?module.exports=w:"function"===typeof define&&define.amd?define("eve",[],function(){return w}):N.eve=w})(this);
(function(N,k){"function"===typeof define&&define.amd?define("Snap.svg",["eve"],function(L){return k(N,L)}):k(N,N.eve)})(this,function(N,k){var L=function(a){var k={},y=N.requestAnimationFrame||N.webkitRequestAnimationFrame||N.mozRequestAnimationFrame||N.oRequestAnimationFrame||N.msRequestAnimationFrame||function(a){setTimeout(a,16)},M=Array.isArray||function(a){return a instanceof Array||"[object Array]"==Object.prototype.toString.call(a)},A=0,w="M"+(+new Date).toString(36),z=function(a){if(null==
a)return this.s;var b=this.s-a;this.b+=this.dur*b;this.B+=this.dur*b;this.s=a},d=function(a){if(null==a)return this.spd;this.spd=a},f=function(a){if(null==a)return this.dur;this.s=this.s*a/this.dur;this.dur=a},n=function(){delete k[this.id];this.update();a("mina.stop."+this.id,this)},u=function(){this.pdif||(delete k[this.id],this.update(),this.pdif=this.get()-this.b)},p=function(){this.pdif&&(this.b=this.get()-this.pdif,delete this.pdif,k[this.id]=this)},b=function(){var a;if(M(this.start)){a=[];
for(var b=0,e=this.start.length;b<e;b++)a[b]=+this.start[b]+(this.end[b]-this.start[b])*this.easing(this.s)}else a=+this.start+(this.end-this.start)*this.easing(this.s);this.set(a)},q=function(){var l=0,b;for(b in k)if(k.hasOwnProperty(b)){var e=k[b],f=e.get();l++;e.s=(f-e.b)/(e.dur/e.spd);1<=e.s&&(delete k[b],e.s=1,l--,function(b){setTimeout(function(){a("mina.finish."+b.id,b)})}(e));e.update()}l&&y(q)},e=function(a,r,s,x,G,h,J){a={id:w+(A++).toString(36),start:a,end:r,b:s,s:0,dur:x-s,spd:1,get:G,
set:h,easing:J||e.linear,status:z,speed:d,duration:f,stop:n,pause:u,resume:p,update:b};k[a.id]=a;r=0;for(var K in k)if(k.hasOwnProperty(K)&&(r++,2==r))break;1==r&&y(q);return a};e.time=Date.now||function(){return+new Date};e.getById=function(a){return k[a]||null};e.linear=function(a){return a};e.easeout=function(a){return Math.pow(a,1.7)};e.easein=function(a){return Math.pow(a,0.48)};e.easeinout=function(a){if(1==a)return 1;if(0==a)return 0;var b=0.48-a/1.04,e=Math.sqrt(0.1734+b*b);a=e-b;a=Math.pow(Math.abs(a),
1/3)*(0>a?-1:1);b=-e-b;b=Math.pow(Math.abs(b),1/3)*(0>b?-1:1);a=a+b+0.5;return 3*(1-a)*a*a+a*a*a};e.backin=function(a){return 1==a?1:a*a*(2.70158*a-1.70158)};e.backout=function(a){if(0==a)return 0;a-=1;return a*a*(2.70158*a+1.70158)+1};e.elastic=function(a){return a==!!a?a:Math.pow(2,-10*a)*Math.sin(2*(a-0.075)*Math.PI/0.3)+1};e.bounce=function(a){a<1/2.75?a*=7.5625*a:a<2/2.75?(a-=1.5/2.75,a=7.5625*a*a+0.75):a<2.5/2.75?(a-=2.25/2.75,a=7.5625*a*a+0.9375):(a-=2.625/2.75,a=7.5625*a*a+0.984375);return a};
return N.mina=e}("undefined"==typeof k?function(){}:k),C=function(){function a(c,t){if(c){if(c.tagName)return x(c);if(y(c,"array")&&a.set)return a.set.apply(a,c);if(c instanceof e)return c;if(null==t)return c=G.doc.querySelector(c),x(c)}return new s(null==c?"100%":c,null==t?"100%":t)}function v(c,a){if(a){"#text"==c&&(c=G.doc.createTextNode(a.text||""));"string"==typeof c&&(c=v(c));if("string"==typeof a)return"xlink:"==a.substring(0,6)?c.getAttributeNS(m,a.substring(6)):"xml:"==a.substring(0,4)?c.getAttributeNS(la,
a.substring(4)):c.getAttribute(a);for(var da in a)if(a[h](da)){var b=J(a[da]);b?"xlink:"==da.substring(0,6)?c.setAttributeNS(m,da.substring(6),b):"xml:"==da.substring(0,4)?c.setAttributeNS(la,da.substring(4),b):c.setAttribute(da,b):c.removeAttribute(da)}}else c=G.doc.createElementNS(la,c);return c}function y(c,a){a=J.prototype.toLowerCase.call(a);return"finite"==a?isFinite(c):"array"==a&&(c instanceof Array||Array.isArray&&Array.isArray(c))?!0:"null"==a&&null===c||a==typeof c&&null!==c||"object"==
a&&c===Object(c)||$.call(c).slice(8,-1).toLowerCase()==a}function M(c){if("function"==typeof c||Object(c)!==c)return c;var a=new c.constructor,b;for(b in c)c[h](b)&&(a[b]=M(c[b]));return a}function A(c,a,b){function m(){var e=Array.prototype.slice.call(arguments,0),f=e.join("\u2400"),d=m.cache=m.cache||{},l=m.count=m.count||[];if(d[h](f)){a:for(var e=l,l=f,B=0,H=e.length;B<H;B++)if(e[B]===l){e.push(e.splice(B,1)[0]);break a}return b?b(d[f]):d[f]}1E3<=l.length&&delete d[l.shift()];l.push(f);d[f]=c.apply(a,
e);return b?b(d[f]):d[f]}return m}function w(c,a,b,m,e,f){return null==e?(c-=b,a-=m,c||a?(180*I.atan2(-a,-c)/C+540)%360:0):w(c,a,e,f)-w(b,m,e,f)}function z(c){return c%360*C/180}function d(c){var a=[];c=c.replace(/(?:^|\s)(\w+)\(([^)]+)\)/g,function(c,b,m){m=m.split(/\s*,\s*|\s+/);"rotate"==b&&1==m.length&&m.push(0,0);"scale"==b&&(2<m.length?m=m.slice(0,2):2==m.length&&m.push(0,0),1==m.length&&m.push(m[0],0,0));"skewX"==b?a.push(["m",1,0,I.tan(z(m[0])),1,0,0]):"skewY"==b?a.push(["m",1,I.tan(z(m[0])),
0,1,0,0]):a.push([b.charAt(0)].concat(m));return c});return a}function f(c,t){var b=O(c),m=new a.Matrix;if(b)for(var e=0,f=b.length;e<f;e++){var h=b[e],d=h.length,B=J(h[0]).toLowerCase(),H=h[0]!=B,l=H?m.invert():0,E;"t"==B&&2==d?m.translate(h[1],0):"t"==B&&3==d?H?(d=l.x(0,0),B=l.y(0,0),H=l.x(h[1],h[2]),l=l.y(h[1],h[2]),m.translate(H-d,l-B)):m.translate(h[1],h[2]):"r"==B?2==d?(E=E||t,m.rotate(h[1],E.x+E.width/2,E.y+E.height/2)):4==d&&(H?(H=l.x(h[2],h[3]),l=l.y(h[2],h[3]),m.rotate(h[1],H,l)):m.rotate(h[1],
h[2],h[3])):"s"==B?2==d||3==d?(E=E||t,m.scale(h[1],h[d-1],E.x+E.width/2,E.y+E.height/2)):4==d?H?(H=l.x(h[2],h[3]),l=l.y(h[2],h[3]),m.scale(h[1],h[1],H,l)):m.scale(h[1],h[1],h[2],h[3]):5==d&&(H?(H=l.x(h[3],h[4]),l=l.y(h[3],h[4]),m.scale(h[1],h[2],H,l)):m.scale(h[1],h[2],h[3],h[4])):"m"==B&&7==d&&m.add(h[1],h[2],h[3],h[4],h[5],h[6])}return m}function n(c,t){if(null==t){var m=!0;t="linearGradient"==c.type||"radialGradient"==c.type?c.node.getAttribute("gradientTransform"):"pattern"==c.type?c.node.getAttribute("patternTransform"):
c.node.getAttribute("transform");if(!t)return new a.Matrix;t=d(t)}else t=a._.rgTransform.test(t)?J(t).replace(/\.{3}|\u2026/g,c._.transform||aa):d(t),y(t,"array")&&(t=a.path?a.path.toString.call(t):J(t)),c._.transform=t;var b=f(t,c.getBBox(1));if(m)return b;c.matrix=b}function u(c){c=c.node.ownerSVGElement&&x(c.node.ownerSVGElement)||c.node.parentNode&&x(c.node.parentNode)||a.select("svg")||a(0,0);var t=c.select("defs"),t=null==t?!1:t.node;t||(t=r("defs",c.node).node);return t}function p(c){return c.node.ownerSVGElement&&
x(c.node.ownerSVGElement)||a.select("svg")}function b(c,a,m){function b(c){if(null==c)return aa;if(c==+c)return c;v(B,{width:c});try{return B.getBBox().width}catch(a){return 0}}function h(c){if(null==c)return aa;if(c==+c)return c;v(B,{height:c});try{return B.getBBox().height}catch(a){return 0}}function e(b,B){null==a?d[b]=B(c.attr(b)||0):b==a&&(d=B(null==m?c.attr(b)||0:m))}var f=p(c).node,d={},B=f.querySelector(".svg---mgr");B||(B=v("rect"),v(B,{x:-9E9,y:-9E9,width:10,height:10,"class":"svg---mgr",
fill:"none"}),f.appendChild(B));switch(c.type){case "rect":e("rx",b),e("ry",h);case "image":e("width",b),e("height",h);case "text":e("x",b);e("y",h);break;case "circle":e("cx",b);e("cy",h);e("r",b);break;case "ellipse":e("cx",b);e("cy",h);e("rx",b);e("ry",h);break;case "line":e("x1",b);e("x2",b);e("y1",h);e("y2",h);break;case "marker":e("refX",b);e("markerWidth",b);e("refY",h);e("markerHeight",h);break;case "radialGradient":e("fx",b);e("fy",h);break;case "tspan":e("dx",b);e("dy",h);break;default:e(a,
b)}f.removeChild(B);return d}function q(c){y(c,"array")||(c=Array.prototype.slice.call(arguments,0));for(var a=0,b=0,m=this.node;this[a];)delete this[a++];for(a=0;a<c.length;a++)"set"==c[a].type?c[a].forEach(function(c){m.appendChild(c.node)}):m.appendChild(c[a].node);for(var h=m.childNodes,a=0;a<h.length;a++)this[b++]=x(h[a]);return this}function e(c){if(c.snap in E)return E[c.snap];var a=this.id=V(),b;try{b=c.ownerSVGElement}catch(m){}this.node=c;b&&(this.paper=new s(b));this.type=c.tagName;this.anims=
{};this._={transform:[]};c.snap=a;E[a]=this;"g"==this.type&&(this.add=q);if(this.type in{g:1,mask:1,pattern:1})for(var e in s.prototype)s.prototype[h](e)&&(this[e]=s.prototype[e])}function l(c){this.node=c}function r(c,a){var b=v(c);a.appendChild(b);return x(b)}function s(c,a){var b,m,f,d=s.prototype;if(c&&"svg"==c.tagName){if(c.snap in E)return E[c.snap];var l=c.ownerDocument;b=new e(c);m=c.getElementsByTagName("desc")[0];f=c.getElementsByTagName("defs")[0];m||(m=v("desc"),m.appendChild(l.createTextNode("Created with Snap")),
b.node.appendChild(m));f||(f=v("defs"),b.node.appendChild(f));b.defs=f;for(var ca in d)d[h](ca)&&(b[ca]=d[ca]);b.paper=b.root=b}else b=r("svg",G.doc.body),v(b.node,{height:a,version:1.1,width:c,xmlns:la});return b}function x(c){return!c||c instanceof e||c instanceof l?c:c.tagName&&"svg"==c.tagName.toLowerCase()?new s(c):c.tagName&&"object"==c.tagName.toLowerCase()&&"image/svg+xml"==c.type?new s(c.contentDocument.getElementsByTagName("svg")[0]):new e(c)}a.version="0.3.0";a.toString=function(){return"Snap v"+
this.version};a._={};var G={win:N,doc:N.document};a._.glob=G;var h="hasOwnProperty",J=String,K=parseFloat,U=parseInt,I=Math,P=I.max,Q=I.min,Y=I.abs,C=I.PI,aa="",$=Object.prototype.toString,F=/^\s*((#[a-f\d]{6})|(#[a-f\d]{3})|rgba?\(\s*([\d\.]+%?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+%?(?:\s*,\s*[\d\.]+%?)?)\s*\)|hsba?\(\s*([\d\.]+(?:deg|\xb0|%)?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+(?:%?\s*,\s*[\d\.]+)?%?)\s*\)|hsla?\(\s*([\d\.]+(?:deg|\xb0|%)?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+(?:%?\s*,\s*[\d\.]+)?%?)\s*\))\s*$/i;a._.separator=
RegExp("[,\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]+");var S=RegExp("[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*"),X={hs:1,rg:1},W=RegExp("([a-z])[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029,]*((-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*)+)",
"ig"),ma=RegExp("([rstm])[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029,]*((-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*)+)","ig"),Z=RegExp("(-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?)[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*",
"ig"),na=0,ba="S"+(+new Date).toString(36),V=function(){return ba+(na++).toString(36)},m="http://www.w3.org/1999/xlink",la="http://www.w3.org/2000/svg",E={},ca=a.url=function(c){return"url('#"+c+"')"};a._.$=v;a._.id=V;a.format=function(){var c=/\{([^\}]+)\}/g,a=/(?:(?:^|\.)(.+?)(?=\[|\.|$|\()|\[('|")(.+?)\2\])(\(\))?/g,b=function(c,b,m){var h=m;b.replace(a,function(c,a,b,m,t){a=a||m;h&&(a in h&&(h=h[a]),"function"==typeof h&&t&&(h=h()))});return h=(null==h||h==m?c:h)+""};return function(a,m){return J(a).replace(c,
function(c,a){return b(c,a,m)})}}();a._.clone=M;a._.cacher=A;a.rad=z;a.deg=function(c){return 180*c/C%360};a.angle=w;a.is=y;a.snapTo=function(c,a,b){b=y(b,"finite")?b:10;if(y(c,"array"))for(var m=c.length;m--;){if(Y(c[m]-a)<=b)return c[m]}else{c=+c;m=a%c;if(m<b)return a-m;if(m>c-b)return a-m+c}return a};a.getRGB=A(function(c){if(!c||(c=J(c)).indexOf("-")+1)return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka};if("none"==c)return{r:-1,g:-1,b:-1,hex:"none",toString:ka};!X[h](c.toLowerCase().substring(0,
2))&&"#"!=c.charAt()&&(c=T(c));if(!c)return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka};var b,m,e,f,d;if(c=c.match(F)){c[2]&&(e=U(c[2].substring(5),16),m=U(c[2].substring(3,5),16),b=U(c[2].substring(1,3),16));c[3]&&(e=U((d=c[3].charAt(3))+d,16),m=U((d=c[3].charAt(2))+d,16),b=U((d=c[3].charAt(1))+d,16));c[4]&&(d=c[4].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b*=2.55),m=K(d[1]),"%"==d[1].slice(-1)&&(m*=2.55),e=K(d[2]),"%"==d[2].slice(-1)&&(e*=2.55),"rgba"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),
d[3]&&"%"==d[3].slice(-1)&&(f/=100));if(c[5])return d=c[5].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b/=100),m=K(d[1]),"%"==d[1].slice(-1)&&(m/=100),e=K(d[2]),"%"==d[2].slice(-1)&&(e/=100),"deg"!=d[0].slice(-3)&&"\u00b0"!=d[0].slice(-1)||(b/=360),"hsba"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),d[3]&&"%"==d[3].slice(-1)&&(f/=100),a.hsb2rgb(b,m,e,f);if(c[6])return d=c[6].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b/=100),m=K(d[1]),"%"==d[1].slice(-1)&&(m/=100),e=K(d[2]),"%"==d[2].slice(-1)&&(e/=100),
"deg"!=d[0].slice(-3)&&"\u00b0"!=d[0].slice(-1)||(b/=360),"hsla"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),d[3]&&"%"==d[3].slice(-1)&&(f/=100),a.hsl2rgb(b,m,e,f);b=Q(I.round(b),255);m=Q(I.round(m),255);e=Q(I.round(e),255);f=Q(P(f,0),1);c={r:b,g:m,b:e,toString:ka};c.hex="#"+(16777216|e|m<<8|b<<16).toString(16).slice(1);c.opacity=y(f,"finite")?f:1;return c}return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka}},a);a.hsb=A(function(c,b,m){return a.hsb2rgb(c,b,m).hex});a.hsl=A(function(c,b,m){return a.hsl2rgb(c,
b,m).hex});a.rgb=A(function(c,a,b,m){if(y(m,"finite")){var e=I.round;return"rgba("+[e(c),e(a),e(b),+m.toFixed(2)]+")"}return"#"+(16777216|b|a<<8|c<<16).toString(16).slice(1)});var T=function(c){var a=G.doc.getElementsByTagName("head")[0]||G.doc.getElementsByTagName("svg")[0];T=A(function(c){if("red"==c.toLowerCase())return"rgb(255, 0, 0)";a.style.color="rgb(255, 0, 0)";a.style.color=c;c=G.doc.defaultView.getComputedStyle(a,aa).getPropertyValue("color");return"rgb(255, 0, 0)"==c?null:c});return T(c)},
qa=function(){return"hsb("+[this.h,this.s,this.b]+")"},ra=function(){return"hsl("+[this.h,this.s,this.l]+")"},ka=function(){return 1==this.opacity||null==this.opacity?this.hex:"rgba("+[this.r,this.g,this.b,this.opacity]+")"},D=function(c,b,m){null==b&&y(c,"object")&&"r"in c&&"g"in c&&"b"in c&&(m=c.b,b=c.g,c=c.r);null==b&&y(c,string)&&(m=a.getRGB(c),c=m.r,b=m.g,m=m.b);if(1<c||1<b||1<m)c/=255,b/=255,m/=255;return[c,b,m]},oa=function(c,b,m,e){c=I.round(255*c);b=I.round(255*b);m=I.round(255*m);c={r:c,
g:b,b:m,opacity:y(e,"finite")?e:1,hex:a.rgb(c,b,m),toString:ka};y(e,"finite")&&(c.opacity=e);return c};a.color=function(c){var b;y(c,"object")&&"h"in c&&"s"in c&&"b"in c?(b=a.hsb2rgb(c),c.r=b.r,c.g=b.g,c.b=b.b,c.opacity=1,c.hex=b.hex):y(c,"object")&&"h"in c&&"s"in c&&"l"in c?(b=a.hsl2rgb(c),c.r=b.r,c.g=b.g,c.b=b.b,c.opacity=1,c.hex=b.hex):(y(c,"string")&&(c=a.getRGB(c)),y(c,"object")&&"r"in c&&"g"in c&&"b"in c&&!("error"in c)?(b=a.rgb2hsl(c),c.h=b.h,c.s=b.s,c.l=b.l,b=a.rgb2hsb(c),c.v=b.b):(c={hex:"none"},
c.r=c.g=c.b=c.h=c.s=c.v=c.l=-1,c.error=1));c.toString=ka;return c};a.hsb2rgb=function(c,a,b,m){y(c,"object")&&"h"in c&&"s"in c&&"b"in c&&(b=c.b,a=c.s,c=c.h,m=c.o);var e,h,d;c=360*c%360/60;d=b*a;a=d*(1-Y(c%2-1));b=e=h=b-d;c=~~c;b+=[d,a,0,0,a,d][c];e+=[a,d,d,a,0,0][c];h+=[0,0,a,d,d,a][c];return oa(b,e,h,m)};a.hsl2rgb=function(c,a,b,m){y(c,"object")&&"h"in c&&"s"in c&&"l"in c&&(b=c.l,a=c.s,c=c.h);if(1<c||1<a||1<b)c/=360,a/=100,b/=100;var e,h,d;c=360*c%360/60;d=2*a*(0.5>b?b:1-b);a=d*(1-Y(c%2-1));b=e=
h=b-d/2;c=~~c;b+=[d,a,0,0,a,d][c];e+=[a,d,d,a,0,0][c];h+=[0,0,a,d,d,a][c];return oa(b,e,h,m)};a.rgb2hsb=function(c,a,b){b=D(c,a,b);c=b[0];a=b[1];b=b[2];var m,e;m=P(c,a,b);e=m-Q(c,a,b);c=((0==e?0:m==c?(a-b)/e:m==a?(b-c)/e+2:(c-a)/e+4)+360)%6*60/360;return{h:c,s:0==e?0:e/m,b:m,toString:qa}};a.rgb2hsl=function(c,a,b){b=D(c,a,b);c=b[0];a=b[1];b=b[2];var m,e,h;m=P(c,a,b);e=Q(c,a,b);h=m-e;c=((0==h?0:m==c?(a-b)/h:m==a?(b-c)/h+2:(c-a)/h+4)+360)%6*60/360;m=(m+e)/2;return{h:c,s:0==h?0:0.5>m?h/(2*m):h/(2-2*
m),l:m,toString:ra}};a.parsePathString=function(c){if(!c)return null;var b=a.path(c);if(b.arr)return a.path.clone(b.arr);var m={a:7,c:6,o:2,h:1,l:2,m:2,r:4,q:4,s:4,t:2,v:1,u:3,z:0},e=[];y(c,"array")&&y(c[0],"array")&&(e=a.path.clone(c));e.length||J(c).replace(W,function(c,a,b){var h=[];c=a.toLowerCase();b.replace(Z,function(c,a){a&&h.push(+a)});"m"==c&&2<h.length&&(e.push([a].concat(h.splice(0,2))),c="l",a="m"==a?"l":"L");"o"==c&&1==h.length&&e.push([a,h[0] ]);if("r"==c)e.push([a].concat(h));else for(;h.length>=
m[c]&&(e.push([a].concat(h.splice(0,m[c]))),m[c]););});e.toString=a.path.toString;b.arr=a.path.clone(e);return e};var O=a.parseTransformString=function(c){if(!c)return null;var b=[];y(c,"array")&&y(c[0],"array")&&(b=a.path.clone(c));b.length||J(c).replace(ma,function(c,a,m){var e=[];a.toLowerCase();m.replace(Z,function(c,a){a&&e.push(+a)});b.push([a].concat(e))});b.toString=a.path.toString;return b};a._.svgTransform2string=d;a._.rgTransform=RegExp("^[a-z][\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*-?\\.?\\d",
"i");a._.transform2matrix=f;a._unit2px=b;a._.getSomeDefs=u;a._.getSomeSVG=p;a.select=function(c){return x(G.doc.querySelector(c))};a.selectAll=function(c){c=G.doc.querySelectorAll(c);for(var b=(a.set||Array)(),m=0;m<c.length;m++)b.push(x(c[m]));return b};setInterval(function(){for(var c in E)if(E[h](c)){var a=E[c],b=a.node;("svg"!=a.type&&!b.ownerSVGElement||"svg"==a.type&&(!b.parentNode||"ownerSVGElement"in b.parentNode&&!b.ownerSVGElement))&&delete E[c]}},1E4);(function(c){function m(c){function a(c,
b){var m=v(c.node,b);(m=(m=m&&m.match(d))&&m[2])&&"#"==m.charAt()&&(m=m.substring(1))&&(f[m]=(f[m]||[]).concat(function(a){var m={};m[b]=ca(a);v(c.node,m)}))}function b(c){var a=v(c.node,"xlink:href");a&&"#"==a.charAt()&&(a=a.substring(1))&&(f[a]=(f[a]||[]).concat(function(a){c.attr("xlink:href","#"+a)}))}var e=c.selectAll("*"),h,d=/^\s*url\(("|'|)(.*)\1\)\s*$/;c=[];for(var f={},l=0,E=e.length;l<E;l++){h=e[l];a(h,"fill");a(h,"stroke");a(h,"filter");a(h,"mask");a(h,"clip-path");b(h);var t=v(h.node,
"id");t&&(v(h.node,{id:h.id}),c.push({old:t,id:h.id}))}l=0;for(E=c.length;l<E;l++)if(e=f[c[l].old])for(h=0,t=e.length;h<t;h++)e[h](c[l].id)}function e(c,a,b){return function(m){m=m.slice(c,a);1==m.length&&(m=m[0]);return b?b(m):m}}function d(c){return function(){var a=c?"<"+this.type:"",b=this.node.attributes,m=this.node.childNodes;if(c)for(var e=0,h=b.length;e<h;e++)a+=" "+b[e].name+'="'+b[e].value.replace(/"/g,'\\"')+'"';if(m.length){c&&(a+=">");e=0;for(h=m.length;e<h;e++)3==m[e].nodeType?a+=m[e].nodeValue:
1==m[e].nodeType&&(a+=x(m[e]).toString());c&&(a+="</"+this.type+">")}else c&&(a+="/>");return a}}c.attr=function(c,a){if(!c)return this;if(y(c,"string"))if(1<arguments.length){var b={};b[c]=a;c=b}else return k("snap.util.getattr."+c,this).firstDefined();for(var m in c)c[h](m)&&k("snap.util.attr."+m,this,c[m]);return this};c.getBBox=function(c){if(!a.Matrix||!a.path)return this.node.getBBox();var b=this,m=new a.Matrix;if(b.removed)return a._.box();for(;"use"==b.type;)if(c||(m=m.add(b.transform().localMatrix.translate(b.attr("x")||
0,b.attr("y")||0))),b.original)b=b.original;else var e=b.attr("xlink:href"),b=b.original=b.node.ownerDocument.getElementById(e.substring(e.indexOf("#")+1));var e=b._,h=a.path.get[b.type]||a.path.get.deflt;try{if(c)return e.bboxwt=h?a.path.getBBox(b.realPath=h(b)):a._.box(b.node.getBBox()),a._.box(e.bboxwt);b.realPath=h(b);b.matrix=b.transform().localMatrix;e.bbox=a.path.getBBox(a.path.map(b.realPath,m.add(b.matrix)));return a._.box(e.bbox)}catch(d){return a._.box()}};var f=function(){return this.string};
c.transform=function(c){var b=this._;if(null==c){var m=this;c=new a.Matrix(this.node.getCTM());for(var e=n(this),h=[e],d=new a.Matrix,l=e.toTransformString(),b=J(e)==J(this.matrix)?J(b.transform):l;"svg"!=m.type&&(m=m.parent());)h.push(n(m));for(m=h.length;m--;)d.add(h[m]);return{string:b,globalMatrix:c,totalMatrix:d,localMatrix:e,diffMatrix:c.clone().add(e.invert()),global:c.toTransformString(),total:d.toTransformString(),local:l,toString:f}}c instanceof a.Matrix?this.matrix=c:n(this,c);this.node&&
("linearGradient"==this.type||"radialGradient"==this.type?v(this.node,{gradientTransform:this.matrix}):"pattern"==this.type?v(this.node,{patternTransform:this.matrix}):v(this.node,{transform:this.matrix}));return this};c.parent=function(){return x(this.node.parentNode)};c.append=c.add=function(c){if(c){if("set"==c.type){var a=this;c.forEach(function(c){a.add(c)});return this}c=x(c);this.node.appendChild(c.node);c.paper=this.paper}return this};c.appendTo=function(c){c&&(c=x(c),c.append(this));return this};
c.prepend=function(c){if(c){if("set"==c.type){var a=this,b;c.forEach(function(c){b?b.after(c):a.prepend(c);b=c});return this}c=x(c);var m=c.parent();this.node.insertBefore(c.node,this.node.firstChild);this.add&&this.add();c.paper=this.paper;this.parent()&&this.parent().add();m&&m.add()}return this};c.prependTo=function(c){c=x(c);c.prepend(this);return this};c.before=function(c){if("set"==c.type){var a=this;c.forEach(function(c){var b=c.parent();a.node.parentNode.insertBefore(c.node,a.node);b&&b.add()});
this.parent().add();return this}c=x(c);var b=c.parent();this.node.parentNode.insertBefore(c.node,this.node);this.parent()&&this.parent().add();b&&b.add();c.paper=this.paper;return this};c.after=function(c){c=x(c);var a=c.parent();this.node.nextSibling?this.node.parentNode.insertBefore(c.node,this.node.nextSibling):this.node.parentNode.appendChild(c.node);this.parent()&&this.parent().add();a&&a.add();c.paper=this.paper;return this};c.insertBefore=function(c){c=x(c);var a=this.parent();c.node.parentNode.insertBefore(this.node,
c.node);this.paper=c.paper;a&&a.add();c.parent()&&c.parent().add();return this};c.insertAfter=function(c){c=x(c);var a=this.parent();c.node.parentNode.insertBefore(this.node,c.node.nextSibling);this.paper=c.paper;a&&a.add();c.parent()&&c.parent().add();return this};c.remove=function(){var c=this.parent();this.node.parentNode&&this.node.parentNode.removeChild(this.node);delete this.paper;this.removed=!0;c&&c.add();return this};c.select=function(c){return x(this.node.querySelector(c))};c.selectAll=
function(c){c=this.node.querySelectorAll(c);for(var b=(a.set||Array)(),m=0;m<c.length;m++)b.push(x(c[m]));return b};c.asPX=function(c,a){null==a&&(a=this.attr(c));return+b(this,c,a)};c.use=function(){var c,a=this.node.id;a||(a=this.id,v(this.node,{id:a}));c="linearGradient"==this.type||"radialGradient"==this.type||"pattern"==this.type?r(this.type,this.node.parentNode):r("use",this.node.parentNode);v(c.node,{"xlink:href":"#"+a});c.original=this;return c};var l=/\S+/g;c.addClass=function(c){var a=(c||
"").match(l)||[];c=this.node;var b=c.className.baseVal,m=b.match(l)||[],e,h,d;if(a.length){for(e=0;d=a[e++];)h=m.indexOf(d),~h||m.push(d);a=m.join(" ");b!=a&&(c.className.baseVal=a)}return this};c.removeClass=function(c){var a=(c||"").match(l)||[];c=this.node;var b=c.className.baseVal,m=b.match(l)||[],e,h;if(m.length){for(e=0;h=a[e++];)h=m.indexOf(h),~h&&m.splice(h,1);a=m.join(" ");b!=a&&(c.className.baseVal=a)}return this};c.hasClass=function(c){return!!~(this.node.className.baseVal.match(l)||[]).indexOf(c)};
c.toggleClass=function(c,a){if(null!=a)return a?this.addClass(c):this.removeClass(c);var b=(c||"").match(l)||[],m=this.node,e=m.className.baseVal,h=e.match(l)||[],d,f,E;for(d=0;E=b[d++];)f=h.indexOf(E),~f?h.splice(f,1):h.push(E);b=h.join(" ");e!=b&&(m.className.baseVal=b);return this};c.clone=function(){var c=x(this.node.cloneNode(!0));v(c.node,"id")&&v(c.node,{id:c.id});m(c);c.insertAfter(this);return c};c.toDefs=function(){u(this).appendChild(this.node);return this};c.pattern=c.toPattern=function(c,
a,b,m){var e=r("pattern",u(this));null==c&&(c=this.getBBox());y(c,"object")&&"x"in c&&(a=c.y,b=c.width,m=c.height,c=c.x);v(e.node,{x:c,y:a,width:b,height:m,patternUnits:"userSpaceOnUse",id:e.id,viewBox:[c,a,b,m].join(" ")});e.node.appendChild(this.node);return e};c.marker=function(c,a,b,m,e,h){var d=r("marker",u(this));null==c&&(c=this.getBBox());y(c,"object")&&"x"in c&&(a=c.y,b=c.width,m=c.height,e=c.refX||c.cx,h=c.refY||c.cy,c=c.x);v(d.node,{viewBox:[c,a,b,m].join(" "),markerWidth:b,markerHeight:m,
orient:"auto",refX:e||0,refY:h||0,id:d.id});d.node.appendChild(this.node);return d};var E=function(c,a,b,m){"function"!=typeof b||b.length||(m=b,b=L.linear);this.attr=c;this.dur=a;b&&(this.easing=b);m&&(this.callback=m)};a._.Animation=E;a.animation=function(c,a,b,m){return new E(c,a,b,m)};c.inAnim=function(){var c=[],a;for(a in this.anims)this.anims[h](a)&&function(a){c.push({anim:new E(a._attrs,a.dur,a.easing,a._callback),mina:a,curStatus:a.status(),status:function(c){return a.status(c)},stop:function(){a.stop()}})}(this.anims[a]);
return c};a.animate=function(c,a,b,m,e,h){"function"!=typeof e||e.length||(h=e,e=L.linear);var d=L.time();c=L(c,a,d,d+m,L.time,b,e);h&&k.once("mina.finish."+c.id,h);return c};c.stop=function(){for(var c=this.inAnim(),a=0,b=c.length;a<b;a++)c[a].stop();return this};c.animate=function(c,a,b,m){"function"!=typeof b||b.length||(m=b,b=L.linear);c instanceof E&&(m=c.callback,b=c.easing,a=b.dur,c=c.attr);var d=[],f=[],l={},t,ca,n,T=this,q;for(q in c)if(c[h](q)){T.equal?(n=T.equal(q,J(c[q])),t=n.from,ca=
n.to,n=n.f):(t=+T.attr(q),ca=+c[q]);var la=y(t,"array")?t.length:1;l[q]=e(d.length,d.length+la,n);d=d.concat(t);f=f.concat(ca)}t=L.time();var p=L(d,f,t,t+a,L.time,function(c){var a={},b;for(b in l)l[h](b)&&(a[b]=l[b](c));T.attr(a)},b);T.anims[p.id]=p;p._attrs=c;p._callback=m;k("snap.animcreated."+T.id,p);k.once("mina.finish."+p.id,function(){delete T.anims[p.id];m&&m.call(T)});k.once("mina.stop."+p.id,function(){delete T.anims[p.id]});return T};var T={};c.data=function(c,b){var m=T[this.id]=T[this.id]||
{};if(0==arguments.length)return k("snap.data.get."+this.id,this,m,null),m;if(1==arguments.length){if(a.is(c,"object")){for(var e in c)c[h](e)&&this.data(e,c[e]);return this}k("snap.data.get."+this.id,this,m[c],c);return m[c]}m[c]=b;k("snap.data.set."+this.id,this,b,c);return this};c.removeData=function(c){null==c?T[this.id]={}:T[this.id]&&delete T[this.id][c];return this};c.outerSVG=c.toString=d(1);c.innerSVG=d()})(e.prototype);a.parse=function(c){var a=G.doc.createDocumentFragment(),b=!0,m=G.doc.createElement("div");
c=J(c);c.match(/^\s*<\s*svg(?:\s|>)/)||(c="<svg>"+c+"</svg>",b=!1);m.innerHTML=c;if(c=m.getElementsByTagName("svg")[0])if(b)a=c;else for(;c.firstChild;)a.appendChild(c.firstChild);m.innerHTML=aa;return new l(a)};l.prototype.select=e.prototype.select;l.prototype.selectAll=e.prototype.selectAll;a.fragment=function(){for(var c=Array.prototype.slice.call(arguments,0),b=G.doc.createDocumentFragment(),m=0,e=c.length;m<e;m++){var h=c[m];h.node&&h.node.nodeType&&b.appendChild(h.node);h.nodeType&&b.appendChild(h);
"string"==typeof h&&b.appendChild(a.parse(h).node)}return new l(b)};a._.make=r;a._.wrap=x;s.prototype.el=function(c,a){var b=r(c,this.node);a&&b.attr(a);return b};k.on("snap.util.getattr",function(){var c=k.nt(),c=c.substring(c.lastIndexOf(".")+1),a=c.replace(/[A-Z]/g,function(c){return"-"+c.toLowerCase()});return pa[h](a)?this.node.ownerDocument.defaultView.getComputedStyle(this.node,null).getPropertyValue(a):v(this.node,c)});var pa={"alignment-baseline":0,"baseline-shift":0,clip:0,"clip-path":0,
"clip-rule":0,color:0,"color-interpolation":0,"color-interpolation-filters":0,"color-profile":0,"color-rendering":0,cursor:0,direction:0,display:0,"dominant-baseline":0,"enable-background":0,fill:0,"fill-opacity":0,"fill-rule":0,filter:0,"flood-color":0,"flood-opacity":0,font:0,"font-family":0,"font-size":0,"font-size-adjust":0,"font-stretch":0,"font-style":0,"font-variant":0,"font-weight":0,"glyph-orientation-horizontal":0,"glyph-orientation-vertical":0,"image-rendering":0,kerning:0,"letter-spacing":0,
"lighting-color":0,marker:0,"marker-end":0,"marker-mid":0,"marker-start":0,mask:0,opacity:0,overflow:0,"pointer-events":0,"shape-rendering":0,"stop-color":0,"stop-opacity":0,stroke:0,"stroke-dasharray":0,"stroke-dashoffset":0,"stroke-linecap":0,"stroke-linejoin":0,"stroke-miterlimit":0,"stroke-opacity":0,"stroke-width":0,"text-anchor":0,"text-decoration":0,"text-rendering":0,"unicode-bidi":0,visibility:0,"word-spacing":0,"writing-mode":0};k.on("snap.util.attr",function(c){var a=k.nt(),b={},a=a.substring(a.lastIndexOf(".")+
1);b[a]=c;var m=a.replace(/-(\w)/gi,function(c,a){return a.toUpperCase()}),a=a.replace(/[A-Z]/g,function(c){return"-"+c.toLowerCase()});pa[h](a)?this.node.style[m]=null==c?aa:c:v(this.node,b)});a.ajax=function(c,a,b,m){var e=new XMLHttpRequest,h=V();if(e){if(y(a,"function"))m=b,b=a,a=null;else if(y(a,"object")){var d=[],f;for(f in a)a.hasOwnProperty(f)&&d.push(encodeURIComponent(f)+"="+encodeURIComponent(a[f]));a=d.join("&")}e.open(a?"POST":"GET",c,!0);a&&(e.setRequestHeader("X-Requested-With","XMLHttpRequest"),
e.setRequestHeader("Content-type","application/x-www-form-urlencoded"));b&&(k.once("snap.ajax."+h+".0",b),k.once("snap.ajax."+h+".200",b),k.once("snap.ajax."+h+".304",b));e.onreadystatechange=function(){4==e.readyState&&k("snap.ajax."+h+"."+e.status,m,e)};if(4==e.readyState)return e;e.send(a);return e}};a.load=function(c,b,m){a.ajax(c,function(c){c=a.parse(c.responseText);m?b.call(m,c):b(c)})};a.getElementByPoint=function(c,a){var b,m,e=G.doc.elementFromPoint(c,a);if(G.win.opera&&"svg"==e.tagName){b=
e;m=b.getBoundingClientRect();b=b.ownerDocument;var h=b.body,d=b.documentElement;b=m.top+(g.win.pageYOffset||d.scrollTop||h.scrollTop)-(d.clientTop||h.clientTop||0);m=m.left+(g.win.pageXOffset||d.scrollLeft||h.scrollLeft)-(d.clientLeft||h.clientLeft||0);h=e.createSVGRect();h.x=c-m;h.y=a-b;h.width=h.height=1;b=e.getIntersectionList(h,null);b.length&&(e=b[b.length-1])}return e?x(e):null};a.plugin=function(c){c(a,e,s,G,l)};return G.win.Snap=a}();C.plugin(function(a,k,y,M,A){function w(a,d,f,b,q,e){null==
d&&"[object SVGMatrix]"==z.call(a)?(this.a=a.a,this.b=a.b,this.c=a.c,this.d=a.d,this.e=a.e,this.f=a.f):null!=a?(this.a=+a,this.b=+d,this.c=+f,this.d=+b,this.e=+q,this.f=+e):(this.a=1,this.c=this.b=0,this.d=1,this.f=this.e=0)}var z=Object.prototype.toString,d=String,f=Math;(function(n){function k(a){return a[0]*a[0]+a[1]*a[1]}function p(a){var d=f.sqrt(k(a));a[0]&&(a[0]/=d);a[1]&&(a[1]/=d)}n.add=function(a,d,e,f,n,p){var k=[[],[],[] ],u=[[this.a,this.c,this.e],[this.b,this.d,this.f],[0,0,1] ];d=[[a,
e,n],[d,f,p],[0,0,1] ];a&&a instanceof w&&(d=[[a.a,a.c,a.e],[a.b,a.d,a.f],[0,0,1] ]);for(a=0;3>a;a++)for(e=0;3>e;e++){for(f=n=0;3>f;f++)n+=u[a][f]*d[f][e];k[a][e]=n}this.a=k[0][0];this.b=k[1][0];this.c=k[0][1];this.d=k[1][1];this.e=k[0][2];this.f=k[1][2];return this};n.invert=function(){var a=this.a*this.d-this.b*this.c;return new w(this.d/a,-this.b/a,-this.c/a,this.a/a,(this.c*this.f-this.d*this.e)/a,(this.b*this.e-this.a*this.f)/a)};n.clone=function(){return new w(this.a,this.b,this.c,this.d,this.e,
this.f)};n.translate=function(a,d){return this.add(1,0,0,1,a,d)};n.scale=function(a,d,e,f){null==d&&(d=a);(e||f)&&this.add(1,0,0,1,e,f);this.add(a,0,0,d,0,0);(e||f)&&this.add(1,0,0,1,-e,-f);return this};n.rotate=function(b,d,e){b=a.rad(b);d=d||0;e=e||0;var l=+f.cos(b).toFixed(9);b=+f.sin(b).toFixed(9);this.add(l,b,-b,l,d,e);return this.add(1,0,0,1,-d,-e)};n.x=function(a,d){return a*this.a+d*this.c+this.e};n.y=function(a,d){return a*this.b+d*this.d+this.f};n.get=function(a){return+this[d.fromCharCode(97+
a)].toFixed(4)};n.toString=function(){return"matrix("+[this.get(0),this.get(1),this.get(2),this.get(3),this.get(4),this.get(5)].join()+")"};n.offset=function(){return[this.e.toFixed(4),this.f.toFixed(4)]};n.determinant=function(){return this.a*this.d-this.b*this.c};n.split=function(){var b={};b.dx=this.e;b.dy=this.f;var d=[[this.a,this.c],[this.b,this.d] ];b.scalex=f.sqrt(k(d[0]));p(d[0]);b.shear=d[0][0]*d[1][0]+d[0][1]*d[1][1];d[1]=[d[1][0]-d[0][0]*b.shear,d[1][1]-d[0][1]*b.shear];b.scaley=f.sqrt(k(d[1]));
p(d[1]);b.shear/=b.scaley;0>this.determinant()&&(b.scalex=-b.scalex);var e=-d[0][1],d=d[1][1];0>d?(b.rotate=a.deg(f.acos(d)),0>e&&(b.rotate=360-b.rotate)):b.rotate=a.deg(f.asin(e));b.isSimple=!+b.shear.toFixed(9)&&(b.scalex.toFixed(9)==b.scaley.toFixed(9)||!b.rotate);b.isSuperSimple=!+b.shear.toFixed(9)&&b.scalex.toFixed(9)==b.scaley.toFixed(9)&&!b.rotate;b.noRotation=!+b.shear.toFixed(9)&&!b.rotate;return b};n.toTransformString=function(a){a=a||this.split();if(+a.shear.toFixed(9))return"m"+[this.get(0),
this.get(1),this.get(2),this.get(3),this.get(4),this.get(5)];a.scalex=+a.scalex.toFixed(4);a.scaley=+a.scaley.toFixed(4);a.rotate=+a.rotate.toFixed(4);return(a.dx||a.dy?"t"+[+a.dx.toFixed(4),+a.dy.toFixed(4)]:"")+(1!=a.scalex||1!=a.scaley?"s"+[a.scalex,a.scaley,0,0]:"")+(a.rotate?"r"+[+a.rotate.toFixed(4),0,0]:"")}})(w.prototype);a.Matrix=w;a.matrix=function(a,d,f,b,k,e){return new w(a,d,f,b,k,e)}});C.plugin(function(a,v,y,M,A){function w(h){return function(d){k.stop();d instanceof A&&1==d.node.childNodes.length&&
("radialGradient"==d.node.firstChild.tagName||"linearGradient"==d.node.firstChild.tagName||"pattern"==d.node.firstChild.tagName)&&(d=d.node.firstChild,b(this).appendChild(d),d=u(d));if(d instanceof v)if("radialGradient"==d.type||"linearGradient"==d.type||"pattern"==d.type){d.node.id||e(d.node,{id:d.id});var f=l(d.node.id)}else f=d.attr(h);else f=a.color(d),f.error?(f=a(b(this).ownerSVGElement).gradient(d))?(f.node.id||e(f.node,{id:f.id}),f=l(f.node.id)):f=d:f=r(f);d={};d[h]=f;e(this.node,d);this.node.style[h]=
x}}function z(a){k.stop();a==+a&&(a+="px");this.node.style.fontSize=a}function d(a){var b=[];a=a.childNodes;for(var e=0,f=a.length;e<f;e++){var l=a[e];3==l.nodeType&&b.push(l.nodeValue);"tspan"==l.tagName&&(1==l.childNodes.length&&3==l.firstChild.nodeType?b.push(l.firstChild.nodeValue):b.push(d(l)))}return b}function f(){k.stop();return this.node.style.fontSize}var n=a._.make,u=a._.wrap,p=a.is,b=a._.getSomeDefs,q=/^url\(#?([^)]+)\)$/,e=a._.$,l=a.url,r=String,s=a._.separator,x="";k.on("snap.util.attr.mask",
function(a){if(a instanceof v||a instanceof A){k.stop();a instanceof A&&1==a.node.childNodes.length&&(a=a.node.firstChild,b(this).appendChild(a),a=u(a));if("mask"==a.type)var d=a;else d=n("mask",b(this)),d.node.appendChild(a.node);!d.node.id&&e(d.node,{id:d.id});e(this.node,{mask:l(d.id)})}});(function(a){k.on("snap.util.attr.clip",a);k.on("snap.util.attr.clip-path",a);k.on("snap.util.attr.clipPath",a)})(function(a){if(a instanceof v||a instanceof A){k.stop();if("clipPath"==a.type)var d=a;else d=
n("clipPath",b(this)),d.node.appendChild(a.node),!d.node.id&&e(d.node,{id:d.id});e(this.node,{"clip-path":l(d.id)})}});k.on("snap.util.attr.fill",w("fill"));k.on("snap.util.attr.stroke",w("stroke"));var G=/^([lr])(?:\(([^)]*)\))?(.*)$/i;k.on("snap.util.grad.parse",function(a){a=r(a);var b=a.match(G);if(!b)return null;a=b[1];var e=b[2],b=b[3],e=e.split(/\s*,\s*/).map(function(a){return+a==a?+a:a});1==e.length&&0==e[0]&&(e=[]);b=b.split("-");b=b.map(function(a){a=a.split(":");var b={color:a[0]};a[1]&&
(b.offset=parseFloat(a[1]));return b});return{type:a,params:e,stops:b}});k.on("snap.util.attr.d",function(b){k.stop();p(b,"array")&&p(b[0],"array")&&(b=a.path.toString.call(b));b=r(b);b.match(/[ruo]/i)&&(b=a.path.toAbsolute(b));e(this.node,{d:b})})(-1);k.on("snap.util.attr.#text",function(a){k.stop();a=r(a);for(a=M.doc.createTextNode(a);this.node.firstChild;)this.node.removeChild(this.node.firstChild);this.node.appendChild(a)})(-1);k.on("snap.util.attr.path",function(a){k.stop();this.attr({d:a})})(-1);
k.on("snap.util.attr.class",function(a){k.stop();this.node.className.baseVal=a})(-1);k.on("snap.util.attr.viewBox",function(a){a=p(a,"object")&&"x"in a?[a.x,a.y,a.width,a.height].join(" "):p(a,"array")?a.join(" "):a;e(this.node,{viewBox:a});k.stop()})(-1);k.on("snap.util.attr.transform",function(a){this.transform(a);k.stop()})(-1);k.on("snap.util.attr.r",function(a){"rect"==this.type&&(k.stop(),e(this.node,{rx:a,ry:a}))})(-1);k.on("snap.util.attr.textpath",function(a){k.stop();if("text"==this.type){var d,
f;if(!a&&this.textPath){for(a=this.textPath;a.node.firstChild;)this.node.appendChild(a.node.firstChild);a.remove();delete this.textPath}else if(p(a,"string")?(d=b(this),a=u(d.parentNode).path(a),d.appendChild(a.node),d=a.id,a.attr({id:d})):(a=u(a),a instanceof v&&(d=a.attr("id"),d||(d=a.id,a.attr({id:d})))),d)if(a=this.textPath,f=this.node,a)a.attr({"xlink:href":"#"+d});else{for(a=e("textPath",{"xlink:href":"#"+d});f.firstChild;)a.appendChild(f.firstChild);f.appendChild(a);this.textPath=u(a)}}})(-1);
k.on("snap.util.attr.text",function(a){if("text"==this.type){for(var b=this.node,d=function(a){var b=e("tspan");if(p(a,"array"))for(var f=0;f<a.length;f++)b.appendChild(d(a[f]));else b.appendChild(M.doc.createTextNode(a));b.normalize&&b.normalize();return b};b.firstChild;)b.removeChild(b.firstChild);for(a=d(a);a.firstChild;)b.appendChild(a.firstChild)}k.stop()})(-1);k.on("snap.util.attr.fontSize",z)(-1);k.on("snap.util.attr.font-size",z)(-1);k.on("snap.util.getattr.transform",function(){k.stop();
return this.transform()})(-1);k.on("snap.util.getattr.textpath",function(){k.stop();return this.textPath})(-1);(function(){function b(d){return function(){k.stop();var b=M.doc.defaultView.getComputedStyle(this.node,null).getPropertyValue("marker-"+d);return"none"==b?b:a(M.doc.getElementById(b.match(q)[1]))}}function d(a){return function(b){k.stop();var d="marker"+a.charAt(0).toUpperCase()+a.substring(1);if(""==b||!b)this.node.style[d]="none";else if("marker"==b.type){var f=b.node.id;f||e(b.node,{id:b.id});
this.node.style[d]=l(f)}}}k.on("snap.util.getattr.marker-end",b("end"))(-1);k.on("snap.util.getattr.markerEnd",b("end"))(-1);k.on("snap.util.getattr.marker-start",b("start"))(-1);k.on("snap.util.getattr.markerStart",b("start"))(-1);k.on("snap.util.getattr.marker-mid",b("mid"))(-1);k.on("snap.util.getattr.markerMid",b("mid"))(-1);k.on("snap.util.attr.marker-end",d("end"))(-1);k.on("snap.util.attr.markerEnd",d("end"))(-1);k.on("snap.util.attr.marker-start",d("start"))(-1);k.on("snap.util.attr.markerStart",
d("start"))(-1);k.on("snap.util.attr.marker-mid",d("mid"))(-1);k.on("snap.util.attr.markerMid",d("mid"))(-1)})();k.on("snap.util.getattr.r",function(){if("rect"==this.type&&e(this.node,"rx")==e(this.node,"ry"))return k.stop(),e(this.node,"rx")})(-1);k.on("snap.util.getattr.text",function(){if("text"==this.type||"tspan"==this.type){k.stop();var a=d(this.node);return 1==a.length?a[0]:a}})(-1);k.on("snap.util.getattr.#text",function(){return this.node.textContent})(-1);k.on("snap.util.getattr.viewBox",
function(){k.stop();var b=e(this.node,"viewBox");if(b)return b=b.split(s),a._.box(+b[0],+b[1],+b[2],+b[3])})(-1);k.on("snap.util.getattr.points",function(){var a=e(this.node,"points");k.stop();if(a)return a.split(s)})(-1);k.on("snap.util.getattr.path",function(){var a=e(this.node,"d");k.stop();return a})(-1);k.on("snap.util.getattr.class",function(){return this.node.className.baseVal})(-1);k.on("snap.util.getattr.fontSize",f)(-1);k.on("snap.util.getattr.font-size",f)(-1)});C.plugin(function(a,v,y,
M,A){function w(a){return a}function z(a){return function(b){return+b.toFixed(3)+a}}var d={"+":function(a,b){return a+b},"-":function(a,b){return a-b},"/":function(a,b){return a/b},"*":function(a,b){return a*b}},f=String,n=/[a-z]+$/i,u=/^\s*([+\-\/*])\s*=\s*([\d.eE+\-]+)\s*([^\d\s]+)?\s*$/;k.on("snap.util.attr",function(a){if(a=f(a).match(u)){var b=k.nt(),b=b.substring(b.lastIndexOf(".")+1),q=this.attr(b),e={};k.stop();var l=a[3]||"",r=q.match(n),s=d[a[1] ];r&&r==l?a=s(parseFloat(q),+a[2]):(q=this.asPX(b),
a=s(this.asPX(b),this.asPX(b,a[2]+l)));isNaN(q)||isNaN(a)||(e[b]=a,this.attr(e))}})(-10);k.on("snap.util.equal",function(a,b){var q=f(this.attr(a)||""),e=f(b).match(u);if(e){k.stop();var l=e[3]||"",r=q.match(n),s=d[e[1] ];if(r&&r==l)return{from:parseFloat(q),to:s(parseFloat(q),+e[2]),f:z(r)};q=this.asPX(a);return{from:q,to:s(q,this.asPX(a,e[2]+l)),f:w}}})(-10)});C.plugin(function(a,v,y,M,A){var w=y.prototype,z=a.is;w.rect=function(a,d,k,p,b,q){var e;null==q&&(q=b);z(a,"object")&&"[object Object]"==
a?e=a:null!=a&&(e={x:a,y:d,width:k,height:p},null!=b&&(e.rx=b,e.ry=q));return this.el("rect",e)};w.circle=function(a,d,k){var p;z(a,"object")&&"[object Object]"==a?p=a:null!=a&&(p={cx:a,cy:d,r:k});return this.el("circle",p)};var d=function(){function a(){this.parentNode.removeChild(this)}return function(d,k){var p=M.doc.createElement("img"),b=M.doc.body;p.style.cssText="position:absolute;left:-9999em;top:-9999em";p.onload=function(){k.call(p);p.onload=p.onerror=null;b.removeChild(p)};p.onerror=a;
b.appendChild(p);p.src=d}}();w.image=function(f,n,k,p,b){var q=this.el("image");if(z(f,"object")&&"src"in f)q.attr(f);else if(null!=f){var e={"xlink:href":f,preserveAspectRatio:"none"};null!=n&&null!=k&&(e.x=n,e.y=k);null!=p&&null!=b?(e.width=p,e.height=b):d(f,function(){a._.$(q.node,{width:this.offsetWidth,height:this.offsetHeight})});a._.$(q.node,e)}return q};w.ellipse=function(a,d,k,p){var b;z(a,"object")&&"[object Object]"==a?b=a:null!=a&&(b={cx:a,cy:d,rx:k,ry:p});return this.el("ellipse",b)};
w.path=function(a){var d;z(a,"object")&&!z(a,"array")?d=a:a&&(d={d:a});return this.el("path",d)};w.group=w.g=function(a){var d=this.el("g");1==arguments.length&&a&&!a.type?d.attr(a):arguments.length&&d.add(Array.prototype.slice.call(arguments,0));return d};w.svg=function(a,d,k,p,b,q,e,l){var r={};z(a,"object")&&null==d?r=a:(null!=a&&(r.x=a),null!=d&&(r.y=d),null!=k&&(r.width=k),null!=p&&(r.height=p),null!=b&&null!=q&&null!=e&&null!=l&&(r.viewBox=[b,q,e,l]));return this.el("svg",r)};w.mask=function(a){var d=
this.el("mask");1==arguments.length&&a&&!a.type?d.attr(a):arguments.length&&d.add(Array.prototype.slice.call(arguments,0));return d};w.ptrn=function(a,d,k,p,b,q,e,l){if(z(a,"object"))var r=a;else arguments.length?(r={},null!=a&&(r.x=a),null!=d&&(r.y=d),null!=k&&(r.width=k),null!=p&&(r.height=p),null!=b&&null!=q&&null!=e&&null!=l&&(r.viewBox=[b,q,e,l])):r={patternUnits:"userSpaceOnUse"};return this.el("pattern",r)};w.use=function(a){return null!=a?(make("use",this.node),a instanceof v&&(a.attr("id")||
a.attr({id:ID()}),a=a.attr("id")),this.el("use",{"xlink:href":a})):v.prototype.use.call(this)};w.text=function(a,d,k){var p={};z(a,"object")?p=a:null!=a&&(p={x:a,y:d,text:k||""});return this.el("text",p)};w.line=function(a,d,k,p){var b={};z(a,"object")?b=a:null!=a&&(b={x1:a,x2:k,y1:d,y2:p});return this.el("line",b)};w.polyline=function(a){1<arguments.length&&(a=Array.prototype.slice.call(arguments,0));var d={};z(a,"object")&&!z(a,"array")?d=a:null!=a&&(d={points:a});return this.el("polyline",d)};
w.polygon=function(a){1<arguments.length&&(a=Array.prototype.slice.call(arguments,0));var d={};z(a,"object")&&!z(a,"array")?d=a:null!=a&&(d={points:a});return this.el("polygon",d)};(function(){function d(){return this.selectAll("stop")}function n(b,d){var f=e("stop"),k={offset:+d+"%"};b=a.color(b);k["stop-color"]=b.hex;1>b.opacity&&(k["stop-opacity"]=b.opacity);e(f,k);this.node.appendChild(f);return this}function u(){if("linearGradient"==this.type){var b=e(this.node,"x1")||0,d=e(this.node,"x2")||
1,f=e(this.node,"y1")||0,k=e(this.node,"y2")||0;return a._.box(b,f,math.abs(d-b),math.abs(k-f))}b=this.node.r||0;return a._.box((this.node.cx||0.5)-b,(this.node.cy||0.5)-b,2*b,2*b)}function p(a,d){function f(a,b){for(var d=(b-u)/(a-w),e=w;e<a;e++)h[e].offset=+(+u+d*(e-w)).toFixed(2);w=a;u=b}var n=k("snap.util.grad.parse",null,d).firstDefined(),p;if(!n)return null;n.params.unshift(a);p="l"==n.type.toLowerCase()?b.apply(0,n.params):q.apply(0,n.params);n.type!=n.type.toLowerCase()&&e(p.node,{gradientUnits:"userSpaceOnUse"});
var h=n.stops,n=h.length,u=0,w=0;n--;for(var v=0;v<n;v++)"offset"in h[v]&&f(v,h[v].offset);h[n].offset=h[n].offset||100;f(n,h[n].offset);for(v=0;v<=n;v++){var y=h[v];p.addStop(y.color,y.offset)}return p}function b(b,k,p,q,w){b=a._.make("linearGradient",b);b.stops=d;b.addStop=n;b.getBBox=u;null!=k&&e(b.node,{x1:k,y1:p,x2:q,y2:w});return b}function q(b,k,p,q,w,h){b=a._.make("radialGradient",b);b.stops=d;b.addStop=n;b.getBBox=u;null!=k&&e(b.node,{cx:k,cy:p,r:q});null!=w&&null!=h&&e(b.node,{fx:w,fy:h});
return b}var e=a._.$;w.gradient=function(a){return p(this.defs,a)};w.gradientLinear=function(a,d,e,f){return b(this.defs,a,d,e,f)};w.gradientRadial=function(a,b,d,e,f){return q(this.defs,a,b,d,e,f)};w.toString=function(){var b=this.node.ownerDocument,d=b.createDocumentFragment(),b=b.createElement("div"),e=this.node.cloneNode(!0);d.appendChild(b);b.appendChild(e);a._.$(e,{xmlns:"http://www.w3.org/2000/svg"});b=b.innerHTML;d.removeChild(d.firstChild);return b};w.clear=function(){for(var a=this.node.firstChild,
b;a;)b=a.nextSibling,"defs"!=a.tagName?a.parentNode.removeChild(a):w.clear.call({node:a}),a=b}})()});C.plugin(function(a,k,y,M){function A(a){var b=A.ps=A.ps||{};b[a]?b[a].sleep=100:b[a]={sleep:100};setTimeout(function(){for(var d in b)b[L](d)&&d!=a&&(b[d].sleep--,!b[d].sleep&&delete b[d])});return b[a]}function w(a,b,d,e){null==a&&(a=b=d=e=0);null==b&&(b=a.y,d=a.width,e=a.height,a=a.x);return{x:a,y:b,width:d,w:d,height:e,h:e,x2:a+d,y2:b+e,cx:a+d/2,cy:b+e/2,r1:F.min(d,e)/2,r2:F.max(d,e)/2,r0:F.sqrt(d*
d+e*e)/2,path:s(a,b,d,e),vb:[a,b,d,e].join(" ")}}function z(){return this.join(",").replace(N,"$1")}function d(a){a=C(a);a.toString=z;return a}function f(a,b,d,h,f,k,l,n,p){if(null==p)return e(a,b,d,h,f,k,l,n);if(0>p||e(a,b,d,h,f,k,l,n)<p)p=void 0;else{var q=0.5,O=1-q,s;for(s=e(a,b,d,h,f,k,l,n,O);0.01<Z(s-p);)q/=2,O+=(s<p?1:-1)*q,s=e(a,b,d,h,f,k,l,n,O);p=O}return u(a,b,d,h,f,k,l,n,p)}function n(b,d){function e(a){return+(+a).toFixed(3)}return a._.cacher(function(a,h,l){a instanceof k&&(a=a.attr("d"));
a=I(a);for(var n,p,D,q,O="",s={},c=0,t=0,r=a.length;t<r;t++){D=a[t];if("M"==D[0])n=+D[1],p=+D[2];else{q=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6]);if(c+q>h){if(d&&!s.start){n=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6],h-c);O+=["C"+e(n.start.x),e(n.start.y),e(n.m.x),e(n.m.y),e(n.x),e(n.y)];if(l)return O;s.start=O;O=["M"+e(n.x),e(n.y)+"C"+e(n.n.x),e(n.n.y),e(n.end.x),e(n.end.y),e(D[5]),e(D[6])].join();c+=q;n=+D[5];p=+D[6];continue}if(!b&&!d)return n=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6],h-c)}c+=q;n=+D[5];p=+D[6]}O+=
D.shift()+D}s.end=O;return n=b?c:d?s:u(n,p,D[0],D[1],D[2],D[3],D[4],D[5],1)},null,a._.clone)}function u(a,b,d,e,h,f,k,l,n){var p=1-n,q=ma(p,3),s=ma(p,2),c=n*n,t=c*n,r=q*a+3*s*n*d+3*p*n*n*h+t*k,q=q*b+3*s*n*e+3*p*n*n*f+t*l,s=a+2*n*(d-a)+c*(h-2*d+a),t=b+2*n*(e-b)+c*(f-2*e+b),x=d+2*n*(h-d)+c*(k-2*h+d),c=e+2*n*(f-e)+c*(l-2*f+e);a=p*a+n*d;b=p*b+n*e;h=p*h+n*k;f=p*f+n*l;l=90-180*F.atan2(s-x,t-c)/S;return{x:r,y:q,m:{x:s,y:t},n:{x:x,y:c},start:{x:a,y:b},end:{x:h,y:f},alpha:l}}function p(b,d,e,h,f,n,k,l){a.is(b,
"array")||(b=[b,d,e,h,f,n,k,l]);b=U.apply(null,b);return w(b.min.x,b.min.y,b.max.x-b.min.x,b.max.y-b.min.y)}function b(a,b,d){return b>=a.x&&b<=a.x+a.width&&d>=a.y&&d<=a.y+a.height}function q(a,d){a=w(a);d=w(d);return b(d,a.x,a.y)||b(d,a.x2,a.y)||b(d,a.x,a.y2)||b(d,a.x2,a.y2)||b(a,d.x,d.y)||b(a,d.x2,d.y)||b(a,d.x,d.y2)||b(a,d.x2,d.y2)||(a.x<d.x2&&a.x>d.x||d.x<a.x2&&d.x>a.x)&&(a.y<d.y2&&a.y>d.y||d.y<a.y2&&d.y>a.y)}function e(a,b,d,e,h,f,n,k,l){null==l&&(l=1);l=(1<l?1:0>l?0:l)/2;for(var p=[-0.1252,
0.1252,-0.3678,0.3678,-0.5873,0.5873,-0.7699,0.7699,-0.9041,0.9041,-0.9816,0.9816],q=[0.2491,0.2491,0.2335,0.2335,0.2032,0.2032,0.1601,0.1601,0.1069,0.1069,0.0472,0.0472],s=0,c=0;12>c;c++)var t=l*p[c]+l,r=t*(t*(-3*a+9*d-9*h+3*n)+6*a-12*d+6*h)-3*a+3*d,t=t*(t*(-3*b+9*e-9*f+3*k)+6*b-12*e+6*f)-3*b+3*e,s=s+q[c]*F.sqrt(r*r+t*t);return l*s}function l(a,b,d){a=I(a);b=I(b);for(var h,f,l,n,k,s,r,O,x,c,t=d?0:[],w=0,v=a.length;w<v;w++)if(x=a[w],"M"==x[0])h=k=x[1],f=s=x[2];else{"C"==x[0]?(x=[h,f].concat(x.slice(1)),
h=x[6],f=x[7]):(x=[h,f,h,f,k,s,k,s],h=k,f=s);for(var G=0,y=b.length;G<y;G++)if(c=b[G],"M"==c[0])l=r=c[1],n=O=c[2];else{"C"==c[0]?(c=[l,n].concat(c.slice(1)),l=c[6],n=c[7]):(c=[l,n,l,n,r,O,r,O],l=r,n=O);var z;var K=x,B=c;z=d;var H=p(K),J=p(B);if(q(H,J)){for(var H=e.apply(0,K),J=e.apply(0,B),H=~~(H/8),J=~~(J/8),U=[],A=[],F={},M=z?0:[],P=0;P<H+1;P++){var C=u.apply(0,K.concat(P/H));U.push({x:C.x,y:C.y,t:P/H})}for(P=0;P<J+1;P++)C=u.apply(0,B.concat(P/J)),A.push({x:C.x,y:C.y,t:P/J});for(P=0;P<H;P++)for(K=
0;K<J;K++){var Q=U[P],L=U[P+1],B=A[K],C=A[K+1],N=0.001>Z(L.x-Q.x)?"y":"x",S=0.001>Z(C.x-B.x)?"y":"x",R;R=Q.x;var Y=Q.y,V=L.x,ea=L.y,fa=B.x,ga=B.y,ha=C.x,ia=C.y;if(W(R,V)<X(fa,ha)||X(R,V)>W(fa,ha)||W(Y,ea)<X(ga,ia)||X(Y,ea)>W(ga,ia))R=void 0;else{var $=(R*ea-Y*V)*(fa-ha)-(R-V)*(fa*ia-ga*ha),aa=(R*ea-Y*V)*(ga-ia)-(Y-ea)*(fa*ia-ga*ha),ja=(R-V)*(ga-ia)-(Y-ea)*(fa-ha);if(ja){var $=$/ja,aa=aa/ja,ja=+$.toFixed(2),ba=+aa.toFixed(2);R=ja<+X(R,V).toFixed(2)||ja>+W(R,V).toFixed(2)||ja<+X(fa,ha).toFixed(2)||
ja>+W(fa,ha).toFixed(2)||ba<+X(Y,ea).toFixed(2)||ba>+W(Y,ea).toFixed(2)||ba<+X(ga,ia).toFixed(2)||ba>+W(ga,ia).toFixed(2)?void 0:{x:$,y:aa}}else R=void 0}R&&F[R.x.toFixed(4)]!=R.y.toFixed(4)&&(F[R.x.toFixed(4)]=R.y.toFixed(4),Q=Q.t+Z((R[N]-Q[N])/(L[N]-Q[N]))*(L.t-Q.t),B=B.t+Z((R[S]-B[S])/(C[S]-B[S]))*(C.t-B.t),0<=Q&&1>=Q&&0<=B&&1>=B&&(z?M++:M.push({x:R.x,y:R.y,t1:Q,t2:B})))}z=M}else z=z?0:[];if(d)t+=z;else{H=0;for(J=z.length;H<J;H++)z[H].segment1=w,z[H].segment2=G,z[H].bez1=x,z[H].bez2=c;t=t.concat(z)}}}return t}
function r(a){var b=A(a);if(b.bbox)return C(b.bbox);if(!a)return w();a=I(a);for(var d=0,e=0,h=[],f=[],l,n=0,k=a.length;n<k;n++)l=a[n],"M"==l[0]?(d=l[1],e=l[2],h.push(d),f.push(e)):(d=U(d,e,l[1],l[2],l[3],l[4],l[5],l[6]),h=h.concat(d.min.x,d.max.x),f=f.concat(d.min.y,d.max.y),d=l[5],e=l[6]);a=X.apply(0,h);l=X.apply(0,f);h=W.apply(0,h);f=W.apply(0,f);f=w(a,l,h-a,f-l);b.bbox=C(f);return f}function s(a,b,d,e,h){if(h)return[["M",+a+ +h,b],["l",d-2*h,0],["a",h,h,0,0,1,h,h],["l",0,e-2*h],["a",h,h,0,0,1,
-h,h],["l",2*h-d,0],["a",h,h,0,0,1,-h,-h],["l",0,2*h-e],["a",h,h,0,0,1,h,-h],["z"] ];a=[["M",a,b],["l",d,0],["l",0,e],["l",-d,0],["z"] ];a.toString=z;return a}function x(a,b,d,e,h){null==h&&null==e&&(e=d);a=+a;b=+b;d=+d;e=+e;if(null!=h){var f=Math.PI/180,l=a+d*Math.cos(-e*f);a+=d*Math.cos(-h*f);var n=b+d*Math.sin(-e*f);b+=d*Math.sin(-h*f);d=[["M",l,n],["A",d,d,0,+(180<h-e),0,a,b] ]}else d=[["M",a,b],["m",0,-e],["a",d,e,0,1,1,0,2*e],["a",d,e,0,1,1,0,-2*e],["z"] ];d.toString=z;return d}function G(b){var e=
A(b);if(e.abs)return d(e.abs);Q(b,"array")&&Q(b&&b[0],"array")||(b=a.parsePathString(b));if(!b||!b.length)return[["M",0,0] ];var h=[],f=0,l=0,n=0,k=0,p=0;"M"==b[0][0]&&(f=+b[0][1],l=+b[0][2],n=f,k=l,p++,h[0]=["M",f,l]);for(var q=3==b.length&&"M"==b[0][0]&&"R"==b[1][0].toUpperCase()&&"Z"==b[2][0].toUpperCase(),s,r,w=p,c=b.length;w<c;w++){h.push(s=[]);r=b[w];p=r[0];if(p!=p.toUpperCase())switch(s[0]=p.toUpperCase(),s[0]){case "A":s[1]=r[1];s[2]=r[2];s[3]=r[3];s[4]=r[4];s[5]=r[5];s[6]=+r[6]+f;s[7]=+r[7]+
l;break;case "V":s[1]=+r[1]+l;break;case "H":s[1]=+r[1]+f;break;case "R":for(var t=[f,l].concat(r.slice(1)),u=2,v=t.length;u<v;u++)t[u]=+t[u]+f,t[++u]=+t[u]+l;h.pop();h=h.concat(P(t,q));break;case "O":h.pop();t=x(f,l,r[1],r[2]);t.push(t[0]);h=h.concat(t);break;case "U":h.pop();h=h.concat(x(f,l,r[1],r[2],r[3]));s=["U"].concat(h[h.length-1].slice(-2));break;case "M":n=+r[1]+f,k=+r[2]+l;default:for(u=1,v=r.length;u<v;u++)s[u]=+r[u]+(u%2?f:l)}else if("R"==p)t=[f,l].concat(r.slice(1)),h.pop(),h=h.concat(P(t,
q)),s=["R"].concat(r.slice(-2));else if("O"==p)h.pop(),t=x(f,l,r[1],r[2]),t.push(t[0]),h=h.concat(t);else if("U"==p)h.pop(),h=h.concat(x(f,l,r[1],r[2],r[3])),s=["U"].concat(h[h.length-1].slice(-2));else for(t=0,u=r.length;t<u;t++)s[t]=r[t];p=p.toUpperCase();if("O"!=p)switch(s[0]){case "Z":f=+n;l=+k;break;case "H":f=s[1];break;case "V":l=s[1];break;case "M":n=s[s.length-2],k=s[s.length-1];default:f=s[s.length-2],l=s[s.length-1]}}h.toString=z;e.abs=d(h);return h}function h(a,b,d,e){return[a,b,d,e,d,
e]}function J(a,b,d,e,h,f){var l=1/3,n=2/3;return[l*a+n*d,l*b+n*e,l*h+n*d,l*f+n*e,h,f]}function K(b,d,e,h,f,l,n,k,p,s){var r=120*S/180,q=S/180*(+f||0),c=[],t,x=a._.cacher(function(a,b,c){var d=a*F.cos(c)-b*F.sin(c);a=a*F.sin(c)+b*F.cos(c);return{x:d,y:a}});if(s)v=s[0],t=s[1],l=s[2],u=s[3];else{t=x(b,d,-q);b=t.x;d=t.y;t=x(k,p,-q);k=t.x;p=t.y;F.cos(S/180*f);F.sin(S/180*f);t=(b-k)/2;v=(d-p)/2;u=t*t/(e*e)+v*v/(h*h);1<u&&(u=F.sqrt(u),e*=u,h*=u);var u=e*e,w=h*h,u=(l==n?-1:1)*F.sqrt(Z((u*w-u*v*v-w*t*t)/
(u*v*v+w*t*t)));l=u*e*v/h+(b+k)/2;var u=u*-h*t/e+(d+p)/2,v=F.asin(((d-u)/h).toFixed(9));t=F.asin(((p-u)/h).toFixed(9));v=b<l?S-v:v;t=k<l?S-t:t;0>v&&(v=2*S+v);0>t&&(t=2*S+t);n&&v>t&&(v-=2*S);!n&&t>v&&(t-=2*S)}if(Z(t-v)>r){var c=t,w=k,G=p;t=v+r*(n&&t>v?1:-1);k=l+e*F.cos(t);p=u+h*F.sin(t);c=K(k,p,e,h,f,0,n,w,G,[t,c,l,u])}l=t-v;f=F.cos(v);r=F.sin(v);n=F.cos(t);t=F.sin(t);l=F.tan(l/4);e=4/3*e*l;l*=4/3*h;h=[b,d];b=[b+e*r,d-l*f];d=[k+e*t,p-l*n];k=[k,p];b[0]=2*h[0]-b[0];b[1]=2*h[1]-b[1];if(s)return[b,d,k].concat(c);
c=[b,d,k].concat(c).join().split(",");s=[];k=0;for(p=c.length;k<p;k++)s[k]=k%2?x(c[k-1],c[k],q).y:x(c[k],c[k+1],q).x;return s}function U(a,b,d,e,h,f,l,k){for(var n=[],p=[[],[] ],s,r,c,t,q=0;2>q;++q)0==q?(r=6*a-12*d+6*h,s=-3*a+9*d-9*h+3*l,c=3*d-3*a):(r=6*b-12*e+6*f,s=-3*b+9*e-9*f+3*k,c=3*e-3*b),1E-12>Z(s)?1E-12>Z(r)||(s=-c/r,0<s&&1>s&&n.push(s)):(t=r*r-4*c*s,c=F.sqrt(t),0>t||(t=(-r+c)/(2*s),0<t&&1>t&&n.push(t),s=(-r-c)/(2*s),0<s&&1>s&&n.push(s)));for(r=q=n.length;q--;)s=n[q],c=1-s,p[0][q]=c*c*c*a+3*
c*c*s*d+3*c*s*s*h+s*s*s*l,p[1][q]=c*c*c*b+3*c*c*s*e+3*c*s*s*f+s*s*s*k;p[0][r]=a;p[1][r]=b;p[0][r+1]=l;p[1][r+1]=k;p[0].length=p[1].length=r+2;return{min:{x:X.apply(0,p[0]),y:X.apply(0,p[1])},max:{x:W.apply(0,p[0]),y:W.apply(0,p[1])}}}function I(a,b){var e=!b&&A(a);if(!b&&e.curve)return d(e.curve);var f=G(a),l=b&&G(b),n={x:0,y:0,bx:0,by:0,X:0,Y:0,qx:null,qy:null},k={x:0,y:0,bx:0,by:0,X:0,Y:0,qx:null,qy:null},p=function(a,b,c){if(!a)return["C",b.x,b.y,b.x,b.y,b.x,b.y];a[0]in{T:1,Q:1}||(b.qx=b.qy=null);
switch(a[0]){case "M":b.X=a[1];b.Y=a[2];break;case "A":a=["C"].concat(K.apply(0,[b.x,b.y].concat(a.slice(1))));break;case "S":"C"==c||"S"==c?(c=2*b.x-b.bx,b=2*b.y-b.by):(c=b.x,b=b.y);a=["C",c,b].concat(a.slice(1));break;case "T":"Q"==c||"T"==c?(b.qx=2*b.x-b.qx,b.qy=2*b.y-b.qy):(b.qx=b.x,b.qy=b.y);a=["C"].concat(J(b.x,b.y,b.qx,b.qy,a[1],a[2]));break;case "Q":b.qx=a[1];b.qy=a[2];a=["C"].concat(J(b.x,b.y,a[1],a[2],a[3],a[4]));break;case "L":a=["C"].concat(h(b.x,b.y,a[1],a[2]));break;case "H":a=["C"].concat(h(b.x,
b.y,a[1],b.y));break;case "V":a=["C"].concat(h(b.x,b.y,b.x,a[1]));break;case "Z":a=["C"].concat(h(b.x,b.y,b.X,b.Y))}return a},s=function(a,b){if(7<a[b].length){a[b].shift();for(var c=a[b];c.length;)q[b]="A",l&&(u[b]="A"),a.splice(b++,0,["C"].concat(c.splice(0,6)));a.splice(b,1);v=W(f.length,l&&l.length||0)}},r=function(a,b,c,d,e){a&&b&&"M"==a[e][0]&&"M"!=b[e][0]&&(b.splice(e,0,["M",d.x,d.y]),c.bx=0,c.by=0,c.x=a[e][1],c.y=a[e][2],v=W(f.length,l&&l.length||0))},q=[],u=[],c="",t="",x=0,v=W(f.length,
l&&l.length||0);for(;x<v;x++){f[x]&&(c=f[x][0]);"C"!=c&&(q[x]=c,x&&(t=q[x-1]));f[x]=p(f[x],n,t);"A"!=q[x]&&"C"==c&&(q[x]="C");s(f,x);l&&(l[x]&&(c=l[x][0]),"C"!=c&&(u[x]=c,x&&(t=u[x-1])),l[x]=p(l[x],k,t),"A"!=u[x]&&"C"==c&&(u[x]="C"),s(l,x));r(f,l,n,k,x);r(l,f,k,n,x);var w=f[x],z=l&&l[x],y=w.length,U=l&&z.length;n.x=w[y-2];n.y=w[y-1];n.bx=$(w[y-4])||n.x;n.by=$(w[y-3])||n.y;k.bx=l&&($(z[U-4])||k.x);k.by=l&&($(z[U-3])||k.y);k.x=l&&z[U-2];k.y=l&&z[U-1]}l||(e.curve=d(f));return l?[f,l]:f}function P(a,
b){for(var d=[],e=0,h=a.length;h-2*!b>e;e+=2){var f=[{x:+a[e-2],y:+a[e-1]},{x:+a[e],y:+a[e+1]},{x:+a[e+2],y:+a[e+3]},{x:+a[e+4],y:+a[e+5]}];b?e?h-4==e?f[3]={x:+a[0],y:+a[1]}:h-2==e&&(f[2]={x:+a[0],y:+a[1]},f[3]={x:+a[2],y:+a[3]}):f[0]={x:+a[h-2],y:+a[h-1]}:h-4==e?f[3]=f[2]:e||(f[0]={x:+a[e],y:+a[e+1]});d.push(["C",(-f[0].x+6*f[1].x+f[2].x)/6,(-f[0].y+6*f[1].y+f[2].y)/6,(f[1].x+6*f[2].x-f[3].x)/6,(f[1].y+6*f[2].y-f[3].y)/6,f[2].x,f[2].y])}return d}y=k.prototype;var Q=a.is,C=a._.clone,L="hasOwnProperty",
N=/,?([a-z]),?/gi,$=parseFloat,F=Math,S=F.PI,X=F.min,W=F.max,ma=F.pow,Z=F.abs;M=n(1);var na=n(),ba=n(0,1),V=a._unit2px;a.path=A;a.path.getTotalLength=M;a.path.getPointAtLength=na;a.path.getSubpath=function(a,b,d){if(1E-6>this.getTotalLength(a)-d)return ba(a,b).end;a=ba(a,d,1);return b?ba(a,b).end:a};y.getTotalLength=function(){if(this.node.getTotalLength)return this.node.getTotalLength()};y.getPointAtLength=function(a){return na(this.attr("d"),a)};y.getSubpath=function(b,d){return a.path.getSubpath(this.attr("d"),
b,d)};a._.box=w;a.path.findDotsAtSegment=u;a.path.bezierBBox=p;a.path.isPointInsideBBox=b;a.path.isBBoxIntersect=q;a.path.intersection=function(a,b){return l(a,b)};a.path.intersectionNumber=function(a,b){return l(a,b,1)};a.path.isPointInside=function(a,d,e){var h=r(a);return b(h,d,e)&&1==l(a,[["M",d,e],["H",h.x2+10] ],1)%2};a.path.getBBox=r;a.path.get={path:function(a){return a.attr("path")},circle:function(a){a=V(a);return x(a.cx,a.cy,a.r)},ellipse:function(a){a=V(a);return x(a.cx||0,a.cy||0,a.rx,
a.ry)},rect:function(a){a=V(a);return s(a.x||0,a.y||0,a.width,a.height,a.rx,a.ry)},image:function(a){a=V(a);return s(a.x||0,a.y||0,a.width,a.height)},line:function(a){return"M"+[a.attr("x1")||0,a.attr("y1")||0,a.attr("x2"),a.attr("y2")]},polyline:function(a){return"M"+a.attr("points")},polygon:function(a){return"M"+a.attr("points")+"z"},deflt:function(a){a=a.node.getBBox();return s(a.x,a.y,a.width,a.height)}};a.path.toRelative=function(b){var e=A(b),h=String.prototype.toLowerCase;if(e.rel)return d(e.rel);
a.is(b,"array")&&a.is(b&&b[0],"array")||(b=a.parsePathString(b));var f=[],l=0,n=0,k=0,p=0,s=0;"M"==b[0][0]&&(l=b[0][1],n=b[0][2],k=l,p=n,s++,f.push(["M",l,n]));for(var r=b.length;s<r;s++){var q=f[s]=[],x=b[s];if(x[0]!=h.call(x[0]))switch(q[0]=h.call(x[0]),q[0]){case "a":q[1]=x[1];q[2]=x[2];q[3]=x[3];q[4]=x[4];q[5]=x[5];q[6]=+(x[6]-l).toFixed(3);q[7]=+(x[7]-n).toFixed(3);break;case "v":q[1]=+(x[1]-n).toFixed(3);break;case "m":k=x[1],p=x[2];default:for(var c=1,t=x.length;c<t;c++)q[c]=+(x[c]-(c%2?l:
n)).toFixed(3)}else for(f[s]=[],"m"==x[0]&&(k=x[1]+l,p=x[2]+n),q=0,c=x.length;q<c;q++)f[s][q]=x[q];x=f[s].length;switch(f[s][0]){case "z":l=k;n=p;break;case "h":l+=+f[s][x-1];break;case "v":n+=+f[s][x-1];break;default:l+=+f[s][x-2],n+=+f[s][x-1]}}f.toString=z;e.rel=d(f);return f};a.path.toAbsolute=G;a.path.toCubic=I;a.path.map=function(a,b){if(!b)return a;var d,e,h,f,l,n,k;a=I(a);h=0;for(l=a.length;h<l;h++)for(k=a[h],f=1,n=k.length;f<n;f+=2)d=b.x(k[f],k[f+1]),e=b.y(k[f],k[f+1]),k[f]=d,k[f+1]=e;return a};
a.path.toString=z;a.path.clone=d});C.plugin(function(a,v,y,C){var A=Math.max,w=Math.min,z=function(a){this.items=[];this.bindings={};this.length=0;this.type="set";if(a)for(var f=0,n=a.length;f<n;f++)a[f]&&(this[this.items.length]=this.items[this.items.length]=a[f],this.length++)};v=z.prototype;v.push=function(){for(var a,f,n=0,k=arguments.length;n<k;n++)if(a=arguments[n])f=this.items.length,this[f]=this.items[f]=a,this.length++;return this};v.pop=function(){this.length&&delete this[this.length--];
return this.items.pop()};v.forEach=function(a,f){for(var n=0,k=this.items.length;n<k&&!1!==a.call(f,this.items[n],n);n++);return this};v.animate=function(d,f,n,u){"function"!=typeof n||n.length||(u=n,n=L.linear);d instanceof a._.Animation&&(u=d.callback,n=d.easing,f=n.dur,d=d.attr);var p=arguments;if(a.is(d,"array")&&a.is(p[p.length-1],"array"))var b=!0;var q,e=function(){q?this.b=q:q=this.b},l=0,r=u&&function(){l++==this.length&&u.call(this)};return this.forEach(function(a,l){k.once("snap.animcreated."+
a.id,e);b?p[l]&&a.animate.apply(a,p[l]):a.animate(d,f,n,r)})};v.remove=function(){for(;this.length;)this.pop().remove();return this};v.bind=function(a,f,k){var u={};if("function"==typeof f)this.bindings[a]=f;else{var p=k||a;this.bindings[a]=function(a){u[p]=a;f.attr(u)}}return this};v.attr=function(a){var f={},k;for(k in a)if(this.bindings[k])this.bindings[k](a[k]);else f[k]=a[k];a=0;for(k=this.items.length;a<k;a++)this.items[a].attr(f);return this};v.clear=function(){for(;this.length;)this.pop()};
v.splice=function(a,f,k){a=0>a?A(this.length+a,0):a;f=A(0,w(this.length-a,f));var u=[],p=[],b=[],q;for(q=2;q<arguments.length;q++)b.push(arguments[q]);for(q=0;q<f;q++)p.push(this[a+q]);for(;q<this.length-a;q++)u.push(this[a+q]);var e=b.length;for(q=0;q<e+u.length;q++)this.items[a+q]=this[a+q]=q<e?b[q]:u[q-e];for(q=this.items.length=this.length-=f-e;this[q];)delete this[q++];return new z(p)};v.exclude=function(a){for(var f=0,k=this.length;f<k;f++)if(this[f]==a)return this.splice(f,1),!0;return!1};
v.insertAfter=function(a){for(var f=this.items.length;f--;)this.items[f].insertAfter(a);return this};v.getBBox=function(){for(var a=[],f=[],k=[],u=[],p=this.items.length;p--;)if(!this.items[p].removed){var b=this.items[p].getBBox();a.push(b.x);f.push(b.y);k.push(b.x+b.width);u.push(b.y+b.height)}a=w.apply(0,a);f=w.apply(0,f);k=A.apply(0,k);u=A.apply(0,u);return{x:a,y:f,x2:k,y2:u,width:k-a,height:u-f,cx:a+(k-a)/2,cy:f+(u-f)/2}};v.clone=function(a){a=new z;for(var f=0,k=this.items.length;f<k;f++)a.push(this.items[f].clone());
return a};v.toString=function(){return"Snap\u2018s set"};v.type="set";a.set=function(){var a=new z;arguments.length&&a.push.apply(a,Array.prototype.slice.call(arguments,0));return a}});C.plugin(function(a,v,y,C){function A(a){var b=a[0];switch(b.toLowerCase()){case "t":return[b,0,0];case "m":return[b,1,0,0,1,0,0];case "r":return 4==a.length?[b,0,a[2],a[3] ]:[b,0];case "s":return 5==a.length?[b,1,1,a[3],a[4] ]:3==a.length?[b,1,1]:[b,1]}}function w(b,d,f){d=q(d).replace(/\.{3}|\u2026/g,b);b=a.parseTransformString(b)||
[];d=a.parseTransformString(d)||[];for(var k=Math.max(b.length,d.length),p=[],v=[],h=0,w,z,y,I;h<k;h++){y=b[h]||A(d[h]);I=d[h]||A(y);if(y[0]!=I[0]||"r"==y[0].toLowerCase()&&(y[2]!=I[2]||y[3]!=I[3])||"s"==y[0].toLowerCase()&&(y[3]!=I[3]||y[4]!=I[4])){b=a._.transform2matrix(b,f());d=a._.transform2matrix(d,f());p=[["m",b.a,b.b,b.c,b.d,b.e,b.f] ];v=[["m",d.a,d.b,d.c,d.d,d.e,d.f] ];break}p[h]=[];v[h]=[];w=0;for(z=Math.max(y.length,I.length);w<z;w++)w in y&&(p[h][w]=y[w]),w in I&&(v[h][w]=I[w])}return{from:u(p),
to:u(v),f:n(p)}}function z(a){return a}function d(a){return function(b){return+b.toFixed(3)+a}}function f(b){return a.rgb(b[0],b[1],b[2])}function n(a){var b=0,d,f,k,n,h,p,q=[];d=0;for(f=a.length;d<f;d++){h="[";p=['"'+a[d][0]+'"'];k=1;for(n=a[d].length;k<n;k++)p[k]="val["+b++ +"]";h+=p+"]";q[d]=h}return Function("val","return Snap.path.toString.call(["+q+"])")}function u(a){for(var b=[],d=0,f=a.length;d<f;d++)for(var k=1,n=a[d].length;k<n;k++)b.push(a[d][k]);return b}var p={},b=/[a-z]+$/i,q=String;
p.stroke=p.fill="colour";v.prototype.equal=function(a,b){return k("snap.util.equal",this,a,b).firstDefined()};k.on("snap.util.equal",function(e,k){var r,s;r=q(this.attr(e)||"");var x=this;if(r==+r&&k==+k)return{from:+r,to:+k,f:z};if("colour"==p[e])return r=a.color(r),s=a.color(k),{from:[r.r,r.g,r.b,r.opacity],to:[s.r,s.g,s.b,s.opacity],f:f};if("transform"==e||"gradientTransform"==e||"patternTransform"==e)return k instanceof a.Matrix&&(k=k.toTransformString()),a._.rgTransform.test(k)||(k=a._.svgTransform2string(k)),
w(r,k,function(){return x.getBBox(1)});if("d"==e||"path"==e)return r=a.path.toCubic(r,k),{from:u(r[0]),to:u(r[1]),f:n(r[0])};if("points"==e)return r=q(r).split(a._.separator),s=q(k).split(a._.separator),{from:r,to:s,f:function(a){return a}};aUnit=r.match(b);s=q(k).match(b);return aUnit&&aUnit==s?{from:parseFloat(r),to:parseFloat(k),f:d(aUnit)}:{from:this.asPX(e),to:this.asPX(e,k),f:z}})});C.plugin(function(a,v,y,C){var A=v.prototype,w="createTouch"in C.doc;v="click dblclick mousedown mousemove mouseout mouseover mouseup touchstart touchmove touchend touchcancel".split(" ");
var z={mousedown:"touchstart",mousemove:"touchmove",mouseup:"touchend"},d=function(a,b){var d="y"==a?"scrollTop":"scrollLeft",e=b&&b.node?b.node.ownerDocument:C.doc;return e[d in e.documentElement?"documentElement":"body"][d]},f=function(){this.returnValue=!1},n=function(){return this.originalEvent.preventDefault()},u=function(){this.cancelBubble=!0},p=function(){return this.originalEvent.stopPropagation()},b=function(){if(C.doc.addEventListener)return function(a,b,e,f){var k=w&&z[b]?z[b]:b,l=function(k){var l=
d("y",f),q=d("x",f);if(w&&z.hasOwnProperty(b))for(var r=0,u=k.targetTouches&&k.targetTouches.length;r<u;r++)if(k.targetTouches[r].target==a||a.contains(k.targetTouches[r].target)){u=k;k=k.targetTouches[r];k.originalEvent=u;k.preventDefault=n;k.stopPropagation=p;break}return e.call(f,k,k.clientX+q,k.clientY+l)};b!==k&&a.addEventListener(b,l,!1);a.addEventListener(k,l,!1);return function(){b!==k&&a.removeEventListener(b,l,!1);a.removeEventListener(k,l,!1);return!0}};if(C.doc.attachEvent)return function(a,
b,e,h){var k=function(a){a=a||h.node.ownerDocument.window.event;var b=d("y",h),k=d("x",h),k=a.clientX+k,b=a.clientY+b;a.preventDefault=a.preventDefault||f;a.stopPropagation=a.stopPropagation||u;return e.call(h,a,k,b)};a.attachEvent("on"+b,k);return function(){a.detachEvent("on"+b,k);return!0}}}(),q=[],e=function(a){for(var b=a.clientX,e=a.clientY,f=d("y"),l=d("x"),n,p=q.length;p--;){n=q[p];if(w)for(var r=a.touches&&a.touches.length,u;r--;){if(u=a.touches[r],u.identifier==n.el._drag.id||n.el.node.contains(u.target)){b=
u.clientX;e=u.clientY;(a.originalEvent?a.originalEvent:a).preventDefault();break}}else a.preventDefault();b+=l;e+=f;k("snap.drag.move."+n.el.id,n.move_scope||n.el,b-n.el._drag.x,e-n.el._drag.y,b,e,a)}},l=function(b){a.unmousemove(e).unmouseup(l);for(var d=q.length,f;d--;)f=q[d],f.el._drag={},k("snap.drag.end."+f.el.id,f.end_scope||f.start_scope||f.move_scope||f.el,b);q=[]};for(y=v.length;y--;)(function(d){a[d]=A[d]=function(e,f){a.is(e,"function")&&(this.events=this.events||[],this.events.push({name:d,
f:e,unbind:b(this.node||document,d,e,f||this)}));return this};a["un"+d]=A["un"+d]=function(a){for(var b=this.events||[],e=b.length;e--;)if(b[e].name==d&&(b[e].f==a||!a)){b[e].unbind();b.splice(e,1);!b.length&&delete this.events;break}return this}})(v[y]);A.hover=function(a,b,d,e){return this.mouseover(a,d).mouseout(b,e||d)};A.unhover=function(a,b){return this.unmouseover(a).unmouseout(b)};var r=[];A.drag=function(b,d,f,h,n,p){function u(r,v,w){(r.originalEvent||r).preventDefault();this._drag.x=v;
this._drag.y=w;this._drag.id=r.identifier;!q.length&&a.mousemove(e).mouseup(l);q.push({el:this,move_scope:h,start_scope:n,end_scope:p});d&&k.on("snap.drag.start."+this.id,d);b&&k.on("snap.drag.move."+this.id,b);f&&k.on("snap.drag.end."+this.id,f);k("snap.drag.start."+this.id,n||h||this,v,w,r)}if(!arguments.length){var v;return this.drag(function(a,b){this.attr({transform:v+(v?"T":"t")+[a,b]})},function(){v=this.transform().local})}this._drag={};r.push({el:this,start:u});this.mousedown(u);return this};
A.undrag=function(){for(var b=r.length;b--;)r[b].el==this&&(this.unmousedown(r[b].start),r.splice(b,1),k.unbind("snap.drag.*."+this.id));!r.length&&a.unmousemove(e).unmouseup(l);return this}});C.plugin(function(a,v,y,C){y=y.prototype;var A=/^\s*url\((.+)\)/,w=String,z=a._.$;a.filter={};y.filter=function(d){var f=this;"svg"!=f.type&&(f=f.paper);d=a.parse(w(d));var k=a._.id(),u=z("filter");z(u,{id:k,filterUnits:"userSpaceOnUse"});u.appendChild(d.node);f.defs.appendChild(u);return new v(u)};k.on("snap.util.getattr.filter",
function(){k.stop();var d=z(this.node,"filter");if(d)return(d=w(d).match(A))&&a.select(d[1])});k.on("snap.util.attr.filter",function(d){if(d instanceof v&&"filter"==d.type){k.stop();var f=d.node.id;f||(z(d.node,{id:d.id}),f=d.id);z(this.node,{filter:a.url(f)})}d&&"none"!=d||(k.stop(),this.node.removeAttribute("filter"))});a.filter.blur=function(d,f){null==d&&(d=2);return a.format('<feGaussianBlur stdDeviation="{def}"/>',{def:null==f?d:[d,f]})};a.filter.blur.toString=function(){return this()};a.filter.shadow=
function(d,f,k,u,p){"string"==typeof k&&(p=u=k,k=4);"string"!=typeof u&&(p=u,u="#000");null==k&&(k=4);null==p&&(p=1);null==d&&(d=0,f=2);null==f&&(f=d);u=a.color(u||"#000");return a.format('<feGaussianBlur in="SourceAlpha" stdDeviation="{blur}"/><feOffset dx="{dx}" dy="{dy}" result="offsetblur"/><feFlood flood-color="{color}"/><feComposite in2="offsetblur" operator="in"/><feComponentTransfer><feFuncA type="linear" slope="{opacity}"/></feComponentTransfer><feMerge><feMergeNode/><feMergeNode in="SourceGraphic"/></feMerge>',
{color:u,dx:d,dy:f,blur:k,opacity:p})};a.filter.shadow.toString=function(){return this()};a.filter.grayscale=function(d){null==d&&(d=1);return a.format('<feColorMatrix type="matrix" values="{a} {b} {c} 0 0 {d} {e} {f} 0 0 {g} {b} {h} 0 0 0 0 0 1 0"/>',{a:0.2126+0.7874*(1-d),b:0.7152-0.7152*(1-d),c:0.0722-0.0722*(1-d),d:0.2126-0.2126*(1-d),e:0.7152+0.2848*(1-d),f:0.0722-0.0722*(1-d),g:0.2126-0.2126*(1-d),h:0.0722+0.9278*(1-d)})};a.filter.grayscale.toString=function(){return this()};a.filter.sepia=
function(d){null==d&&(d=1);return a.format('<feColorMatrix type="matrix" values="{a} {b} {c} 0 0 {d} {e} {f} 0 0 {g} {h} {i} 0 0 0 0 0 1 0"/>',{a:0.393+0.607*(1-d),b:0.769-0.769*(1-d),c:0.189-0.189*(1-d),d:0.349-0.349*(1-d),e:0.686+0.314*(1-d),f:0.168-0.168*(1-d),g:0.272-0.272*(1-d),h:0.534-0.534*(1-d),i:0.131+0.869*(1-d)})};a.filter.sepia.toString=function(){return this()};a.filter.saturate=function(d){null==d&&(d=1);return a.format('<feColorMatrix type="saturate" values="{amount}"/>',{amount:1-
d})};a.filter.saturate.toString=function(){return this()};a.filter.hueRotate=function(d){return a.format('<feColorMatrix type="hueRotate" values="{angle}"/>',{angle:d||0})};a.filter.hueRotate.toString=function(){return this()};a.filter.invert=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="table" tableValues="{amount} {amount2}"/><feFuncG type="table" tableValues="{amount} {amount2}"/><feFuncB type="table" tableValues="{amount} {amount2}"/></feComponentTransfer>',{amount:d,
amount2:1-d})};a.filter.invert.toString=function(){return this()};a.filter.brightness=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="linear" slope="{amount}"/><feFuncG type="linear" slope="{amount}"/><feFuncB type="linear" slope="{amount}"/></feComponentTransfer>',{amount:d})};a.filter.brightness.toString=function(){return this()};a.filter.contrast=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="linear" slope="{amount}" intercept="{amount2}"/><feFuncG type="linear" slope="{amount}" intercept="{amount2}"/><feFuncB type="linear" slope="{amount}" intercept="{amount2}"/></feComponentTransfer>',
{amount:d,amount2:0.5-d/2})};a.filter.contrast.toString=function(){return this()}});return C});

]]> </script>
</svg>



<div class="markdown"><h3>Blocks for d-level systems</h3>
</div>

<pre class='language-julia'><code class='language-julia'>block_3level = matblock(rand_unitary(3); nlevel=3, tag="3 level")</code></pre>
<pre id='var-block_3level' class='code-output documenter-example-output'>3 level</pre>

<pre class='language-julia'><code class='language-julia'>mat(put(4, 2=&gt;block_3level))</code></pre>
<pre id='var-hash107447' class='code-output documenter-example-output'>81×81 SparseMatrixCSC{ComplexF64, Int64} with 243 stored entries:
⢕⢕⢕⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠑⠑⠑⢄⢄⢄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⢕⢕⢕⠅⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠁⠁⢕⢕⢼⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠒⠗⠕⢅⢀⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⢐⢕⢕⠅⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠁⢑⢔⢔⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠰⠵⠕⢅⢀⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢐⢕⢕⠅⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠁⢑⢔⢔⢄⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠐⢕⢕⢕⢀⢀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢐⢕⢕⢕⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⠑⠑⢄⢆⢆
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠨⢕⢕⢕</pre>

<pre class='language-julia'><code class='language-julia'>vizcircuit(put(4, 2=&gt;block_3level))</code></pre>
<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:xlink="http://www.w3.org/1999/xlink"
     version="1.2"
     width="42.5mm" height="75mm" viewBox="0 0 42.5 75"
     stroke="none"
     fill="#000000"
     stroke-width="0.3"
     font-size="3.88"

     id="img-9471da4a">
<defs>
  <marker id="arrow" markerWidth="15" markerHeight="7" refX="5" refY="3.5" orient="auto" markerUnits="strokeWidth">
    <path d="M0,0 L15,3.5 L0,7 z" stroke="context-stroke" fill="context-stroke"/>
  </marker>
</defs>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-9471da4a-1">
  <g transform="translate(23.12,25)">
    <path d="M-8.31,-5 L8.31,-5 8.31,5 -8.31,5  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-9471da4a-2">
  <g transform="translate(23.12,62.5)">
    <path fill="none" d="M10.62,0 L-10.62,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-9471da4a-3">
  <g transform="translate(23.12,43.75)">
    <path fill="none" d="M10.62,0 L-10.62,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-9471da4a-4">
  <g transform="translate(33.75,25)">
    <path fill="none" d="M0,0 L0,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-9471da4a-5">
  <g transform="translate(23.12,6.25)">
    <path fill="none" d="M10.62,0 L-10.62,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-9471da4a-6">
  <g transform="translate(32.59,25)">
    <path fill="none" d="M-1.16,0 L1.16,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-9471da4a-7">
  <g transform="translate(13.66,25)">
    <path fill="none" d="M-1.16,0 L1.16,0 " class="primitive"/>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="3.53" id="img-9471da4a-8">
  <g transform="translate(23.12,25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">3 level</text>
    </g>
  </g>
</g>
<script> <![CDATA[
(function(N){var k=/[\.\/]/,L=/\s*,\s*/,C=function(a,d){return a-d},a,v,y={n:{}},M=function(){for(var a=0,d=this.length;a<d;a++)if("undefined"!=typeof this[a])return this[a]},A=function(){for(var a=this.length;--a;)if("undefined"!=typeof this[a])return this[a]},w=function(k,d){k=String(k);var f=v,n=Array.prototype.slice.call(arguments,2),u=w.listeners(k),p=0,b,q=[],e={},l=[],r=a;l.firstDefined=M;l.lastDefined=A;a=k;for(var s=v=0,x=u.length;s<x;s++)"zIndex"in u[s]&&(q.push(u[s].zIndex),0>u[s].zIndex&&
(e[u[s].zIndex]=u[s]));for(q.sort(C);0>q[p];)if(b=e[q[p++] ],l.push(b.apply(d,n)),v)return v=f,l;for(s=0;s<x;s++)if(b=u[s],"zIndex"in b)if(b.zIndex==q[p]){l.push(b.apply(d,n));if(v)break;do if(p++,(b=e[q[p] ])&&l.push(b.apply(d,n)),v)break;while(b)}else e[b.zIndex]=b;else if(l.push(b.apply(d,n)),v)break;v=f;a=r;return l};w._events=y;w.listeners=function(a){a=a.split(k);var d=y,f,n,u,p,b,q,e,l=[d],r=[];u=0;for(p=a.length;u<p;u++){e=[];b=0;for(q=l.length;b<q;b++)for(d=l[b].n,f=[d[a[u] ],d["*"] ],n=2;n--;)if(d=
f[n])e.push(d),r=r.concat(d.f||[]);l=e}return r};w.on=function(a,d){a=String(a);if("function"!=typeof d)return function(){};for(var f=a.split(L),n=0,u=f.length;n<u;n++)(function(a){a=a.split(k);for(var b=y,f,e=0,l=a.length;e<l;e++)b=b.n,b=b.hasOwnProperty(a[e])&&b[a[e] ]||(b[a[e] ]={n:{}});b.f=b.f||[];e=0;for(l=b.f.length;e<l;e++)if(b.f[e]==d){f=!0;break}!f&&b.f.push(d)})(f[n]);return function(a){+a==+a&&(d.zIndex=+a)}};w.f=function(a){var d=[].slice.call(arguments,1);return function(){w.apply(null,
[a,null].concat(d).concat([].slice.call(arguments,0)))}};w.stop=function(){v=1};w.nt=function(k){return k?(new RegExp("(?:\\.|\\/|^)"+k+"(?:\\.|\\/|$)")).test(a):a};w.nts=function(){return a.split(k)};w.off=w.unbind=function(a,d){if(a){var f=a.split(L);if(1<f.length)for(var n=0,u=f.length;n<u;n++)w.off(f[n],d);else{for(var f=a.split(k),p,b,q,e,l=[y],n=0,u=f.length;n<u;n++)for(e=0;e<l.length;e+=q.length-2){q=[e,1];p=l[e].n;if("*"!=f[n])p[f[n] ]&&q.push(p[f[n] ]);else for(b in p)p.hasOwnProperty(b)&&
q.push(p[b]);l.splice.apply(l,q)}n=0;for(u=l.length;n<u;n++)for(p=l[n];p.n;){if(d){if(p.f){e=0;for(f=p.f.length;e<f;e++)if(p.f[e]==d){p.f.splice(e,1);break}!p.f.length&&delete p.f}for(b in p.n)if(p.n.hasOwnProperty(b)&&p.n[b].f){q=p.n[b].f;e=0;for(f=q.length;e<f;e++)if(q[e]==d){q.splice(e,1);break}!q.length&&delete p.n[b].f}}else for(b in delete p.f,p.n)p.n.hasOwnProperty(b)&&p.n[b].f&&delete p.n[b].f;p=p.n}}}else w._events=y={n:{}}};w.once=function(a,d){var f=function(){w.unbind(a,f);return d.apply(this,
arguments)};return w.on(a,f)};w.version="0.4.2";w.toString=function(){return"You are running Eve 0.4.2"};"undefined"!=typeof module&&module.exports?module.exports=w:"function"===typeof define&&define.amd?define("eve",[],function(){return w}):N.eve=w})(this);
(function(N,k){"function"===typeof define&&define.amd?define("Snap.svg",["eve"],function(L){return k(N,L)}):k(N,N.eve)})(this,function(N,k){var L=function(a){var k={},y=N.requestAnimationFrame||N.webkitRequestAnimationFrame||N.mozRequestAnimationFrame||N.oRequestAnimationFrame||N.msRequestAnimationFrame||function(a){setTimeout(a,16)},M=Array.isArray||function(a){return a instanceof Array||"[object Array]"==Object.prototype.toString.call(a)},A=0,w="M"+(+new Date).toString(36),z=function(a){if(null==
a)return this.s;var b=this.s-a;this.b+=this.dur*b;this.B+=this.dur*b;this.s=a},d=function(a){if(null==a)return this.spd;this.spd=a},f=function(a){if(null==a)return this.dur;this.s=this.s*a/this.dur;this.dur=a},n=function(){delete k[this.id];this.update();a("mina.stop."+this.id,this)},u=function(){this.pdif||(delete k[this.id],this.update(),this.pdif=this.get()-this.b)},p=function(){this.pdif&&(this.b=this.get()-this.pdif,delete this.pdif,k[this.id]=this)},b=function(){var a;if(M(this.start)){a=[];
for(var b=0,e=this.start.length;b<e;b++)a[b]=+this.start[b]+(this.end[b]-this.start[b])*this.easing(this.s)}else a=+this.start+(this.end-this.start)*this.easing(this.s);this.set(a)},q=function(){var l=0,b;for(b in k)if(k.hasOwnProperty(b)){var e=k[b],f=e.get();l++;e.s=(f-e.b)/(e.dur/e.spd);1<=e.s&&(delete k[b],e.s=1,l--,function(b){setTimeout(function(){a("mina.finish."+b.id,b)})}(e));e.update()}l&&y(q)},e=function(a,r,s,x,G,h,J){a={id:w+(A++).toString(36),start:a,end:r,b:s,s:0,dur:x-s,spd:1,get:G,
set:h,easing:J||e.linear,status:z,speed:d,duration:f,stop:n,pause:u,resume:p,update:b};k[a.id]=a;r=0;for(var K in k)if(k.hasOwnProperty(K)&&(r++,2==r))break;1==r&&y(q);return a};e.time=Date.now||function(){return+new Date};e.getById=function(a){return k[a]||null};e.linear=function(a){return a};e.easeout=function(a){return Math.pow(a,1.7)};e.easein=function(a){return Math.pow(a,0.48)};e.easeinout=function(a){if(1==a)return 1;if(0==a)return 0;var b=0.48-a/1.04,e=Math.sqrt(0.1734+b*b);a=e-b;a=Math.pow(Math.abs(a),
1/3)*(0>a?-1:1);b=-e-b;b=Math.pow(Math.abs(b),1/3)*(0>b?-1:1);a=a+b+0.5;return 3*(1-a)*a*a+a*a*a};e.backin=function(a){return 1==a?1:a*a*(2.70158*a-1.70158)};e.backout=function(a){if(0==a)return 0;a-=1;return a*a*(2.70158*a+1.70158)+1};e.elastic=function(a){return a==!!a?a:Math.pow(2,-10*a)*Math.sin(2*(a-0.075)*Math.PI/0.3)+1};e.bounce=function(a){a<1/2.75?a*=7.5625*a:a<2/2.75?(a-=1.5/2.75,a=7.5625*a*a+0.75):a<2.5/2.75?(a-=2.25/2.75,a=7.5625*a*a+0.9375):(a-=2.625/2.75,a=7.5625*a*a+0.984375);return a};
return N.mina=e}("undefined"==typeof k?function(){}:k),C=function(){function a(c,t){if(c){if(c.tagName)return x(c);if(y(c,"array")&&a.set)return a.set.apply(a,c);if(c instanceof e)return c;if(null==t)return c=G.doc.querySelector(c),x(c)}return new s(null==c?"100%":c,null==t?"100%":t)}function v(c,a){if(a){"#text"==c&&(c=G.doc.createTextNode(a.text||""));"string"==typeof c&&(c=v(c));if("string"==typeof a)return"xlink:"==a.substring(0,6)?c.getAttributeNS(m,a.substring(6)):"xml:"==a.substring(0,4)?c.getAttributeNS(la,
a.substring(4)):c.getAttribute(a);for(var da in a)if(a[h](da)){var b=J(a[da]);b?"xlink:"==da.substring(0,6)?c.setAttributeNS(m,da.substring(6),b):"xml:"==da.substring(0,4)?c.setAttributeNS(la,da.substring(4),b):c.setAttribute(da,b):c.removeAttribute(da)}}else c=G.doc.createElementNS(la,c);return c}function y(c,a){a=J.prototype.toLowerCase.call(a);return"finite"==a?isFinite(c):"array"==a&&(c instanceof Array||Array.isArray&&Array.isArray(c))?!0:"null"==a&&null===c||a==typeof c&&null!==c||"object"==
a&&c===Object(c)||$.call(c).slice(8,-1).toLowerCase()==a}function M(c){if("function"==typeof c||Object(c)!==c)return c;var a=new c.constructor,b;for(b in c)c[h](b)&&(a[b]=M(c[b]));return a}function A(c,a,b){function m(){var e=Array.prototype.slice.call(arguments,0),f=e.join("\u2400"),d=m.cache=m.cache||{},l=m.count=m.count||[];if(d[h](f)){a:for(var e=l,l=f,B=0,H=e.length;B<H;B++)if(e[B]===l){e.push(e.splice(B,1)[0]);break a}return b?b(d[f]):d[f]}1E3<=l.length&&delete d[l.shift()];l.push(f);d[f]=c.apply(a,
e);return b?b(d[f]):d[f]}return m}function w(c,a,b,m,e,f){return null==e?(c-=b,a-=m,c||a?(180*I.atan2(-a,-c)/C+540)%360:0):w(c,a,e,f)-w(b,m,e,f)}function z(c){return c%360*C/180}function d(c){var a=[];c=c.replace(/(?:^|\s)(\w+)\(([^)]+)\)/g,function(c,b,m){m=m.split(/\s*,\s*|\s+/);"rotate"==b&&1==m.length&&m.push(0,0);"scale"==b&&(2<m.length?m=m.slice(0,2):2==m.length&&m.push(0,0),1==m.length&&m.push(m[0],0,0));"skewX"==b?a.push(["m",1,0,I.tan(z(m[0])),1,0,0]):"skewY"==b?a.push(["m",1,I.tan(z(m[0])),
0,1,0,0]):a.push([b.charAt(0)].concat(m));return c});return a}function f(c,t){var b=O(c),m=new a.Matrix;if(b)for(var e=0,f=b.length;e<f;e++){var h=b[e],d=h.length,B=J(h[0]).toLowerCase(),H=h[0]!=B,l=H?m.invert():0,E;"t"==B&&2==d?m.translate(h[1],0):"t"==B&&3==d?H?(d=l.x(0,0),B=l.y(0,0),H=l.x(h[1],h[2]),l=l.y(h[1],h[2]),m.translate(H-d,l-B)):m.translate(h[1],h[2]):"r"==B?2==d?(E=E||t,m.rotate(h[1],E.x+E.width/2,E.y+E.height/2)):4==d&&(H?(H=l.x(h[2],h[3]),l=l.y(h[2],h[3]),m.rotate(h[1],H,l)):m.rotate(h[1],
h[2],h[3])):"s"==B?2==d||3==d?(E=E||t,m.scale(h[1],h[d-1],E.x+E.width/2,E.y+E.height/2)):4==d?H?(H=l.x(h[2],h[3]),l=l.y(h[2],h[3]),m.scale(h[1],h[1],H,l)):m.scale(h[1],h[1],h[2],h[3]):5==d&&(H?(H=l.x(h[3],h[4]),l=l.y(h[3],h[4]),m.scale(h[1],h[2],H,l)):m.scale(h[1],h[2],h[3],h[4])):"m"==B&&7==d&&m.add(h[1],h[2],h[3],h[4],h[5],h[6])}return m}function n(c,t){if(null==t){var m=!0;t="linearGradient"==c.type||"radialGradient"==c.type?c.node.getAttribute("gradientTransform"):"pattern"==c.type?c.node.getAttribute("patternTransform"):
c.node.getAttribute("transform");if(!t)return new a.Matrix;t=d(t)}else t=a._.rgTransform.test(t)?J(t).replace(/\.{3}|\u2026/g,c._.transform||aa):d(t),y(t,"array")&&(t=a.path?a.path.toString.call(t):J(t)),c._.transform=t;var b=f(t,c.getBBox(1));if(m)return b;c.matrix=b}function u(c){c=c.node.ownerSVGElement&&x(c.node.ownerSVGElement)||c.node.parentNode&&x(c.node.parentNode)||a.select("svg")||a(0,0);var t=c.select("defs"),t=null==t?!1:t.node;t||(t=r("defs",c.node).node);return t}function p(c){return c.node.ownerSVGElement&&
x(c.node.ownerSVGElement)||a.select("svg")}function b(c,a,m){function b(c){if(null==c)return aa;if(c==+c)return c;v(B,{width:c});try{return B.getBBox().width}catch(a){return 0}}function h(c){if(null==c)return aa;if(c==+c)return c;v(B,{height:c});try{return B.getBBox().height}catch(a){return 0}}function e(b,B){null==a?d[b]=B(c.attr(b)||0):b==a&&(d=B(null==m?c.attr(b)||0:m))}var f=p(c).node,d={},B=f.querySelector(".svg---mgr");B||(B=v("rect"),v(B,{x:-9E9,y:-9E9,width:10,height:10,"class":"svg---mgr",
fill:"none"}),f.appendChild(B));switch(c.type){case "rect":e("rx",b),e("ry",h);case "image":e("width",b),e("height",h);case "text":e("x",b);e("y",h);break;case "circle":e("cx",b);e("cy",h);e("r",b);break;case "ellipse":e("cx",b);e("cy",h);e("rx",b);e("ry",h);break;case "line":e("x1",b);e("x2",b);e("y1",h);e("y2",h);break;case "marker":e("refX",b);e("markerWidth",b);e("refY",h);e("markerHeight",h);break;case "radialGradient":e("fx",b);e("fy",h);break;case "tspan":e("dx",b);e("dy",h);break;default:e(a,
b)}f.removeChild(B);return d}function q(c){y(c,"array")||(c=Array.prototype.slice.call(arguments,0));for(var a=0,b=0,m=this.node;this[a];)delete this[a++];for(a=0;a<c.length;a++)"set"==c[a].type?c[a].forEach(function(c){m.appendChild(c.node)}):m.appendChild(c[a].node);for(var h=m.childNodes,a=0;a<h.length;a++)this[b++]=x(h[a]);return this}function e(c){if(c.snap in E)return E[c.snap];var a=this.id=V(),b;try{b=c.ownerSVGElement}catch(m){}this.node=c;b&&(this.paper=new s(b));this.type=c.tagName;this.anims=
{};this._={transform:[]};c.snap=a;E[a]=this;"g"==this.type&&(this.add=q);if(this.type in{g:1,mask:1,pattern:1})for(var e in s.prototype)s.prototype[h](e)&&(this[e]=s.prototype[e])}function l(c){this.node=c}function r(c,a){var b=v(c);a.appendChild(b);return x(b)}function s(c,a){var b,m,f,d=s.prototype;if(c&&"svg"==c.tagName){if(c.snap in E)return E[c.snap];var l=c.ownerDocument;b=new e(c);m=c.getElementsByTagName("desc")[0];f=c.getElementsByTagName("defs")[0];m||(m=v("desc"),m.appendChild(l.createTextNode("Created with Snap")),
b.node.appendChild(m));f||(f=v("defs"),b.node.appendChild(f));b.defs=f;for(var ca in d)d[h](ca)&&(b[ca]=d[ca]);b.paper=b.root=b}else b=r("svg",G.doc.body),v(b.node,{height:a,version:1.1,width:c,xmlns:la});return b}function x(c){return!c||c instanceof e||c instanceof l?c:c.tagName&&"svg"==c.tagName.toLowerCase()?new s(c):c.tagName&&"object"==c.tagName.toLowerCase()&&"image/svg+xml"==c.type?new s(c.contentDocument.getElementsByTagName("svg")[0]):new e(c)}a.version="0.3.0";a.toString=function(){return"Snap v"+
this.version};a._={};var G={win:N,doc:N.document};a._.glob=G;var h="hasOwnProperty",J=String,K=parseFloat,U=parseInt,I=Math,P=I.max,Q=I.min,Y=I.abs,C=I.PI,aa="",$=Object.prototype.toString,F=/^\s*((#[a-f\d]{6})|(#[a-f\d]{3})|rgba?\(\s*([\d\.]+%?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+%?(?:\s*,\s*[\d\.]+%?)?)\s*\)|hsba?\(\s*([\d\.]+(?:deg|\xb0|%)?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+(?:%?\s*,\s*[\d\.]+)?%?)\s*\)|hsla?\(\s*([\d\.]+(?:deg|\xb0|%)?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+(?:%?\s*,\s*[\d\.]+)?%?)\s*\))\s*$/i;a._.separator=
RegExp("[,\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]+");var S=RegExp("[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*"),X={hs:1,rg:1},W=RegExp("([a-z])[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029,]*((-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*)+)",
"ig"),ma=RegExp("([rstm])[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029,]*((-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*)+)","ig"),Z=RegExp("(-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?)[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*",
"ig"),na=0,ba="S"+(+new Date).toString(36),V=function(){return ba+(na++).toString(36)},m="http://www.w3.org/1999/xlink",la="http://www.w3.org/2000/svg",E={},ca=a.url=function(c){return"url('#"+c+"')"};a._.$=v;a._.id=V;a.format=function(){var c=/\{([^\}]+)\}/g,a=/(?:(?:^|\.)(.+?)(?=\[|\.|$|\()|\[('|")(.+?)\2\])(\(\))?/g,b=function(c,b,m){var h=m;b.replace(a,function(c,a,b,m,t){a=a||m;h&&(a in h&&(h=h[a]),"function"==typeof h&&t&&(h=h()))});return h=(null==h||h==m?c:h)+""};return function(a,m){return J(a).replace(c,
function(c,a){return b(c,a,m)})}}();a._.clone=M;a._.cacher=A;a.rad=z;a.deg=function(c){return 180*c/C%360};a.angle=w;a.is=y;a.snapTo=function(c,a,b){b=y(b,"finite")?b:10;if(y(c,"array"))for(var m=c.length;m--;){if(Y(c[m]-a)<=b)return c[m]}else{c=+c;m=a%c;if(m<b)return a-m;if(m>c-b)return a-m+c}return a};a.getRGB=A(function(c){if(!c||(c=J(c)).indexOf("-")+1)return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka};if("none"==c)return{r:-1,g:-1,b:-1,hex:"none",toString:ka};!X[h](c.toLowerCase().substring(0,
2))&&"#"!=c.charAt()&&(c=T(c));if(!c)return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka};var b,m,e,f,d;if(c=c.match(F)){c[2]&&(e=U(c[2].substring(5),16),m=U(c[2].substring(3,5),16),b=U(c[2].substring(1,3),16));c[3]&&(e=U((d=c[3].charAt(3))+d,16),m=U((d=c[3].charAt(2))+d,16),b=U((d=c[3].charAt(1))+d,16));c[4]&&(d=c[4].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b*=2.55),m=K(d[1]),"%"==d[1].slice(-1)&&(m*=2.55),e=K(d[2]),"%"==d[2].slice(-1)&&(e*=2.55),"rgba"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),
d[3]&&"%"==d[3].slice(-1)&&(f/=100));if(c[5])return d=c[5].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b/=100),m=K(d[1]),"%"==d[1].slice(-1)&&(m/=100),e=K(d[2]),"%"==d[2].slice(-1)&&(e/=100),"deg"!=d[0].slice(-3)&&"\u00b0"!=d[0].slice(-1)||(b/=360),"hsba"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),d[3]&&"%"==d[3].slice(-1)&&(f/=100),a.hsb2rgb(b,m,e,f);if(c[6])return d=c[6].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b/=100),m=K(d[1]),"%"==d[1].slice(-1)&&(m/=100),e=K(d[2]),"%"==d[2].slice(-1)&&(e/=100),
"deg"!=d[0].slice(-3)&&"\u00b0"!=d[0].slice(-1)||(b/=360),"hsla"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),d[3]&&"%"==d[3].slice(-1)&&(f/=100),a.hsl2rgb(b,m,e,f);b=Q(I.round(b),255);m=Q(I.round(m),255);e=Q(I.round(e),255);f=Q(P(f,0),1);c={r:b,g:m,b:e,toString:ka};c.hex="#"+(16777216|e|m<<8|b<<16).toString(16).slice(1);c.opacity=y(f,"finite")?f:1;return c}return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka}},a);a.hsb=A(function(c,b,m){return a.hsb2rgb(c,b,m).hex});a.hsl=A(function(c,b,m){return a.hsl2rgb(c,
b,m).hex});a.rgb=A(function(c,a,b,m){if(y(m,"finite")){var e=I.round;return"rgba("+[e(c),e(a),e(b),+m.toFixed(2)]+")"}return"#"+(16777216|b|a<<8|c<<16).toString(16).slice(1)});var T=function(c){var a=G.doc.getElementsByTagName("head")[0]||G.doc.getElementsByTagName("svg")[0];T=A(function(c){if("red"==c.toLowerCase())return"rgb(255, 0, 0)";a.style.color="rgb(255, 0, 0)";a.style.color=c;c=G.doc.defaultView.getComputedStyle(a,aa).getPropertyValue("color");return"rgb(255, 0, 0)"==c?null:c});return T(c)},
qa=function(){return"hsb("+[this.h,this.s,this.b]+")"},ra=function(){return"hsl("+[this.h,this.s,this.l]+")"},ka=function(){return 1==this.opacity||null==this.opacity?this.hex:"rgba("+[this.r,this.g,this.b,this.opacity]+")"},D=function(c,b,m){null==b&&y(c,"object")&&"r"in c&&"g"in c&&"b"in c&&(m=c.b,b=c.g,c=c.r);null==b&&y(c,string)&&(m=a.getRGB(c),c=m.r,b=m.g,m=m.b);if(1<c||1<b||1<m)c/=255,b/=255,m/=255;return[c,b,m]},oa=function(c,b,m,e){c=I.round(255*c);b=I.round(255*b);m=I.round(255*m);c={r:c,
g:b,b:m,opacity:y(e,"finite")?e:1,hex:a.rgb(c,b,m),toString:ka};y(e,"finite")&&(c.opacity=e);return c};a.color=function(c){var b;y(c,"object")&&"h"in c&&"s"in c&&"b"in c?(b=a.hsb2rgb(c),c.r=b.r,c.g=b.g,c.b=b.b,c.opacity=1,c.hex=b.hex):y(c,"object")&&"h"in c&&"s"in c&&"l"in c?(b=a.hsl2rgb(c),c.r=b.r,c.g=b.g,c.b=b.b,c.opacity=1,c.hex=b.hex):(y(c,"string")&&(c=a.getRGB(c)),y(c,"object")&&"r"in c&&"g"in c&&"b"in c&&!("error"in c)?(b=a.rgb2hsl(c),c.h=b.h,c.s=b.s,c.l=b.l,b=a.rgb2hsb(c),c.v=b.b):(c={hex:"none"},
c.r=c.g=c.b=c.h=c.s=c.v=c.l=-1,c.error=1));c.toString=ka;return c};a.hsb2rgb=function(c,a,b,m){y(c,"object")&&"h"in c&&"s"in c&&"b"in c&&(b=c.b,a=c.s,c=c.h,m=c.o);var e,h,d;c=360*c%360/60;d=b*a;a=d*(1-Y(c%2-1));b=e=h=b-d;c=~~c;b+=[d,a,0,0,a,d][c];e+=[a,d,d,a,0,0][c];h+=[0,0,a,d,d,a][c];return oa(b,e,h,m)};a.hsl2rgb=function(c,a,b,m){y(c,"object")&&"h"in c&&"s"in c&&"l"in c&&(b=c.l,a=c.s,c=c.h);if(1<c||1<a||1<b)c/=360,a/=100,b/=100;var e,h,d;c=360*c%360/60;d=2*a*(0.5>b?b:1-b);a=d*(1-Y(c%2-1));b=e=
h=b-d/2;c=~~c;b+=[d,a,0,0,a,d][c];e+=[a,d,d,a,0,0][c];h+=[0,0,a,d,d,a][c];return oa(b,e,h,m)};a.rgb2hsb=function(c,a,b){b=D(c,a,b);c=b[0];a=b[1];b=b[2];var m,e;m=P(c,a,b);e=m-Q(c,a,b);c=((0==e?0:m==c?(a-b)/e:m==a?(b-c)/e+2:(c-a)/e+4)+360)%6*60/360;return{h:c,s:0==e?0:e/m,b:m,toString:qa}};a.rgb2hsl=function(c,a,b){b=D(c,a,b);c=b[0];a=b[1];b=b[2];var m,e,h;m=P(c,a,b);e=Q(c,a,b);h=m-e;c=((0==h?0:m==c?(a-b)/h:m==a?(b-c)/h+2:(c-a)/h+4)+360)%6*60/360;m=(m+e)/2;return{h:c,s:0==h?0:0.5>m?h/(2*m):h/(2-2*
m),l:m,toString:ra}};a.parsePathString=function(c){if(!c)return null;var b=a.path(c);if(b.arr)return a.path.clone(b.arr);var m={a:7,c:6,o:2,h:1,l:2,m:2,r:4,q:4,s:4,t:2,v:1,u:3,z:0},e=[];y(c,"array")&&y(c[0],"array")&&(e=a.path.clone(c));e.length||J(c).replace(W,function(c,a,b){var h=[];c=a.toLowerCase();b.replace(Z,function(c,a){a&&h.push(+a)});"m"==c&&2<h.length&&(e.push([a].concat(h.splice(0,2))),c="l",a="m"==a?"l":"L");"o"==c&&1==h.length&&e.push([a,h[0] ]);if("r"==c)e.push([a].concat(h));else for(;h.length>=
m[c]&&(e.push([a].concat(h.splice(0,m[c]))),m[c]););});e.toString=a.path.toString;b.arr=a.path.clone(e);return e};var O=a.parseTransformString=function(c){if(!c)return null;var b=[];y(c,"array")&&y(c[0],"array")&&(b=a.path.clone(c));b.length||J(c).replace(ma,function(c,a,m){var e=[];a.toLowerCase();m.replace(Z,function(c,a){a&&e.push(+a)});b.push([a].concat(e))});b.toString=a.path.toString;return b};a._.svgTransform2string=d;a._.rgTransform=RegExp("^[a-z][\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*-?\\.?\\d",
"i");a._.transform2matrix=f;a._unit2px=b;a._.getSomeDefs=u;a._.getSomeSVG=p;a.select=function(c){return x(G.doc.querySelector(c))};a.selectAll=function(c){c=G.doc.querySelectorAll(c);for(var b=(a.set||Array)(),m=0;m<c.length;m++)b.push(x(c[m]));return b};setInterval(function(){for(var c in E)if(E[h](c)){var a=E[c],b=a.node;("svg"!=a.type&&!b.ownerSVGElement||"svg"==a.type&&(!b.parentNode||"ownerSVGElement"in b.parentNode&&!b.ownerSVGElement))&&delete E[c]}},1E4);(function(c){function m(c){function a(c,
b){var m=v(c.node,b);(m=(m=m&&m.match(d))&&m[2])&&"#"==m.charAt()&&(m=m.substring(1))&&(f[m]=(f[m]||[]).concat(function(a){var m={};m[b]=ca(a);v(c.node,m)}))}function b(c){var a=v(c.node,"xlink:href");a&&"#"==a.charAt()&&(a=a.substring(1))&&(f[a]=(f[a]||[]).concat(function(a){c.attr("xlink:href","#"+a)}))}var e=c.selectAll("*"),h,d=/^\s*url\(("|'|)(.*)\1\)\s*$/;c=[];for(var f={},l=0,E=e.length;l<E;l++){h=e[l];a(h,"fill");a(h,"stroke");a(h,"filter");a(h,"mask");a(h,"clip-path");b(h);var t=v(h.node,
"id");t&&(v(h.node,{id:h.id}),c.push({old:t,id:h.id}))}l=0;for(E=c.length;l<E;l++)if(e=f[c[l].old])for(h=0,t=e.length;h<t;h++)e[h](c[l].id)}function e(c,a,b){return function(m){m=m.slice(c,a);1==m.length&&(m=m[0]);return b?b(m):m}}function d(c){return function(){var a=c?"<"+this.type:"",b=this.node.attributes,m=this.node.childNodes;if(c)for(var e=0,h=b.length;e<h;e++)a+=" "+b[e].name+'="'+b[e].value.replace(/"/g,'\\"')+'"';if(m.length){c&&(a+=">");e=0;for(h=m.length;e<h;e++)3==m[e].nodeType?a+=m[e].nodeValue:
1==m[e].nodeType&&(a+=x(m[e]).toString());c&&(a+="</"+this.type+">")}else c&&(a+="/>");return a}}c.attr=function(c,a){if(!c)return this;if(y(c,"string"))if(1<arguments.length){var b={};b[c]=a;c=b}else return k("snap.util.getattr."+c,this).firstDefined();for(var m in c)c[h](m)&&k("snap.util.attr."+m,this,c[m]);return this};c.getBBox=function(c){if(!a.Matrix||!a.path)return this.node.getBBox();var b=this,m=new a.Matrix;if(b.removed)return a._.box();for(;"use"==b.type;)if(c||(m=m.add(b.transform().localMatrix.translate(b.attr("x")||
0,b.attr("y")||0))),b.original)b=b.original;else var e=b.attr("xlink:href"),b=b.original=b.node.ownerDocument.getElementById(e.substring(e.indexOf("#")+1));var e=b._,h=a.path.get[b.type]||a.path.get.deflt;try{if(c)return e.bboxwt=h?a.path.getBBox(b.realPath=h(b)):a._.box(b.node.getBBox()),a._.box(e.bboxwt);b.realPath=h(b);b.matrix=b.transform().localMatrix;e.bbox=a.path.getBBox(a.path.map(b.realPath,m.add(b.matrix)));return a._.box(e.bbox)}catch(d){return a._.box()}};var f=function(){return this.string};
c.transform=function(c){var b=this._;if(null==c){var m=this;c=new a.Matrix(this.node.getCTM());for(var e=n(this),h=[e],d=new a.Matrix,l=e.toTransformString(),b=J(e)==J(this.matrix)?J(b.transform):l;"svg"!=m.type&&(m=m.parent());)h.push(n(m));for(m=h.length;m--;)d.add(h[m]);return{string:b,globalMatrix:c,totalMatrix:d,localMatrix:e,diffMatrix:c.clone().add(e.invert()),global:c.toTransformString(),total:d.toTransformString(),local:l,toString:f}}c instanceof a.Matrix?this.matrix=c:n(this,c);this.node&&
("linearGradient"==this.type||"radialGradient"==this.type?v(this.node,{gradientTransform:this.matrix}):"pattern"==this.type?v(this.node,{patternTransform:this.matrix}):v(this.node,{transform:this.matrix}));return this};c.parent=function(){return x(this.node.parentNode)};c.append=c.add=function(c){if(c){if("set"==c.type){var a=this;c.forEach(function(c){a.add(c)});return this}c=x(c);this.node.appendChild(c.node);c.paper=this.paper}return this};c.appendTo=function(c){c&&(c=x(c),c.append(this));return this};
c.prepend=function(c){if(c){if("set"==c.type){var a=this,b;c.forEach(function(c){b?b.after(c):a.prepend(c);b=c});return this}c=x(c);var m=c.parent();this.node.insertBefore(c.node,this.node.firstChild);this.add&&this.add();c.paper=this.paper;this.parent()&&this.parent().add();m&&m.add()}return this};c.prependTo=function(c){c=x(c);c.prepend(this);return this};c.before=function(c){if("set"==c.type){var a=this;c.forEach(function(c){var b=c.parent();a.node.parentNode.insertBefore(c.node,a.node);b&&b.add()});
this.parent().add();return this}c=x(c);var b=c.parent();this.node.parentNode.insertBefore(c.node,this.node);this.parent()&&this.parent().add();b&&b.add();c.paper=this.paper;return this};c.after=function(c){c=x(c);var a=c.parent();this.node.nextSibling?this.node.parentNode.insertBefore(c.node,this.node.nextSibling):this.node.parentNode.appendChild(c.node);this.parent()&&this.parent().add();a&&a.add();c.paper=this.paper;return this};c.insertBefore=function(c){c=x(c);var a=this.parent();c.node.parentNode.insertBefore(this.node,
c.node);this.paper=c.paper;a&&a.add();c.parent()&&c.parent().add();return this};c.insertAfter=function(c){c=x(c);var a=this.parent();c.node.parentNode.insertBefore(this.node,c.node.nextSibling);this.paper=c.paper;a&&a.add();c.parent()&&c.parent().add();return this};c.remove=function(){var c=this.parent();this.node.parentNode&&this.node.parentNode.removeChild(this.node);delete this.paper;this.removed=!0;c&&c.add();return this};c.select=function(c){return x(this.node.querySelector(c))};c.selectAll=
function(c){c=this.node.querySelectorAll(c);for(var b=(a.set||Array)(),m=0;m<c.length;m++)b.push(x(c[m]));return b};c.asPX=function(c,a){null==a&&(a=this.attr(c));return+b(this,c,a)};c.use=function(){var c,a=this.node.id;a||(a=this.id,v(this.node,{id:a}));c="linearGradient"==this.type||"radialGradient"==this.type||"pattern"==this.type?r(this.type,this.node.parentNode):r("use",this.node.parentNode);v(c.node,{"xlink:href":"#"+a});c.original=this;return c};var l=/\S+/g;c.addClass=function(c){var a=(c||
"").match(l)||[];c=this.node;var b=c.className.baseVal,m=b.match(l)||[],e,h,d;if(a.length){for(e=0;d=a[e++];)h=m.indexOf(d),~h||m.push(d);a=m.join(" ");b!=a&&(c.className.baseVal=a)}return this};c.removeClass=function(c){var a=(c||"").match(l)||[];c=this.node;var b=c.className.baseVal,m=b.match(l)||[],e,h;if(m.length){for(e=0;h=a[e++];)h=m.indexOf(h),~h&&m.splice(h,1);a=m.join(" ");b!=a&&(c.className.baseVal=a)}return this};c.hasClass=function(c){return!!~(this.node.className.baseVal.match(l)||[]).indexOf(c)};
c.toggleClass=function(c,a){if(null!=a)return a?this.addClass(c):this.removeClass(c);var b=(c||"").match(l)||[],m=this.node,e=m.className.baseVal,h=e.match(l)||[],d,f,E;for(d=0;E=b[d++];)f=h.indexOf(E),~f?h.splice(f,1):h.push(E);b=h.join(" ");e!=b&&(m.className.baseVal=b);return this};c.clone=function(){var c=x(this.node.cloneNode(!0));v(c.node,"id")&&v(c.node,{id:c.id});m(c);c.insertAfter(this);return c};c.toDefs=function(){u(this).appendChild(this.node);return this};c.pattern=c.toPattern=function(c,
a,b,m){var e=r("pattern",u(this));null==c&&(c=this.getBBox());y(c,"object")&&"x"in c&&(a=c.y,b=c.width,m=c.height,c=c.x);v(e.node,{x:c,y:a,width:b,height:m,patternUnits:"userSpaceOnUse",id:e.id,viewBox:[c,a,b,m].join(" ")});e.node.appendChild(this.node);return e};c.marker=function(c,a,b,m,e,h){var d=r("marker",u(this));null==c&&(c=this.getBBox());y(c,"object")&&"x"in c&&(a=c.y,b=c.width,m=c.height,e=c.refX||c.cx,h=c.refY||c.cy,c=c.x);v(d.node,{viewBox:[c,a,b,m].join(" "),markerWidth:b,markerHeight:m,
orient:"auto",refX:e||0,refY:h||0,id:d.id});d.node.appendChild(this.node);return d};var E=function(c,a,b,m){"function"!=typeof b||b.length||(m=b,b=L.linear);this.attr=c;this.dur=a;b&&(this.easing=b);m&&(this.callback=m)};a._.Animation=E;a.animation=function(c,a,b,m){return new E(c,a,b,m)};c.inAnim=function(){var c=[],a;for(a in this.anims)this.anims[h](a)&&function(a){c.push({anim:new E(a._attrs,a.dur,a.easing,a._callback),mina:a,curStatus:a.status(),status:function(c){return a.status(c)},stop:function(){a.stop()}})}(this.anims[a]);
return c};a.animate=function(c,a,b,m,e,h){"function"!=typeof e||e.length||(h=e,e=L.linear);var d=L.time();c=L(c,a,d,d+m,L.time,b,e);h&&k.once("mina.finish."+c.id,h);return c};c.stop=function(){for(var c=this.inAnim(),a=0,b=c.length;a<b;a++)c[a].stop();return this};c.animate=function(c,a,b,m){"function"!=typeof b||b.length||(m=b,b=L.linear);c instanceof E&&(m=c.callback,b=c.easing,a=b.dur,c=c.attr);var d=[],f=[],l={},t,ca,n,T=this,q;for(q in c)if(c[h](q)){T.equal?(n=T.equal(q,J(c[q])),t=n.from,ca=
n.to,n=n.f):(t=+T.attr(q),ca=+c[q]);var la=y(t,"array")?t.length:1;l[q]=e(d.length,d.length+la,n);d=d.concat(t);f=f.concat(ca)}t=L.time();var p=L(d,f,t,t+a,L.time,function(c){var a={},b;for(b in l)l[h](b)&&(a[b]=l[b](c));T.attr(a)},b);T.anims[p.id]=p;p._attrs=c;p._callback=m;k("snap.animcreated."+T.id,p);k.once("mina.finish."+p.id,function(){delete T.anims[p.id];m&&m.call(T)});k.once("mina.stop."+p.id,function(){delete T.anims[p.id]});return T};var T={};c.data=function(c,b){var m=T[this.id]=T[this.id]||
{};if(0==arguments.length)return k("snap.data.get."+this.id,this,m,null),m;if(1==arguments.length){if(a.is(c,"object")){for(var e in c)c[h](e)&&this.data(e,c[e]);return this}k("snap.data.get."+this.id,this,m[c],c);return m[c]}m[c]=b;k("snap.data.set."+this.id,this,b,c);return this};c.removeData=function(c){null==c?T[this.id]={}:T[this.id]&&delete T[this.id][c];return this};c.outerSVG=c.toString=d(1);c.innerSVG=d()})(e.prototype);a.parse=function(c){var a=G.doc.createDocumentFragment(),b=!0,m=G.doc.createElement("div");
c=J(c);c.match(/^\s*<\s*svg(?:\s|>)/)||(c="<svg>"+c+"</svg>",b=!1);m.innerHTML=c;if(c=m.getElementsByTagName("svg")[0])if(b)a=c;else for(;c.firstChild;)a.appendChild(c.firstChild);m.innerHTML=aa;return new l(a)};l.prototype.select=e.prototype.select;l.prototype.selectAll=e.prototype.selectAll;a.fragment=function(){for(var c=Array.prototype.slice.call(arguments,0),b=G.doc.createDocumentFragment(),m=0,e=c.length;m<e;m++){var h=c[m];h.node&&h.node.nodeType&&b.appendChild(h.node);h.nodeType&&b.appendChild(h);
"string"==typeof h&&b.appendChild(a.parse(h).node)}return new l(b)};a._.make=r;a._.wrap=x;s.prototype.el=function(c,a){var b=r(c,this.node);a&&b.attr(a);return b};k.on("snap.util.getattr",function(){var c=k.nt(),c=c.substring(c.lastIndexOf(".")+1),a=c.replace(/[A-Z]/g,function(c){return"-"+c.toLowerCase()});return pa[h](a)?this.node.ownerDocument.defaultView.getComputedStyle(this.node,null).getPropertyValue(a):v(this.node,c)});var pa={"alignment-baseline":0,"baseline-shift":0,clip:0,"clip-path":0,
"clip-rule":0,color:0,"color-interpolation":0,"color-interpolation-filters":0,"color-profile":0,"color-rendering":0,cursor:0,direction:0,display:0,"dominant-baseline":0,"enable-background":0,fill:0,"fill-opacity":0,"fill-rule":0,filter:0,"flood-color":0,"flood-opacity":0,font:0,"font-family":0,"font-size":0,"font-size-adjust":0,"font-stretch":0,"font-style":0,"font-variant":0,"font-weight":0,"glyph-orientation-horizontal":0,"glyph-orientation-vertical":0,"image-rendering":0,kerning:0,"letter-spacing":0,
"lighting-color":0,marker:0,"marker-end":0,"marker-mid":0,"marker-start":0,mask:0,opacity:0,overflow:0,"pointer-events":0,"shape-rendering":0,"stop-color":0,"stop-opacity":0,stroke:0,"stroke-dasharray":0,"stroke-dashoffset":0,"stroke-linecap":0,"stroke-linejoin":0,"stroke-miterlimit":0,"stroke-opacity":0,"stroke-width":0,"text-anchor":0,"text-decoration":0,"text-rendering":0,"unicode-bidi":0,visibility:0,"word-spacing":0,"writing-mode":0};k.on("snap.util.attr",function(c){var a=k.nt(),b={},a=a.substring(a.lastIndexOf(".")+
1);b[a]=c;var m=a.replace(/-(\w)/gi,function(c,a){return a.toUpperCase()}),a=a.replace(/[A-Z]/g,function(c){return"-"+c.toLowerCase()});pa[h](a)?this.node.style[m]=null==c?aa:c:v(this.node,b)});a.ajax=function(c,a,b,m){var e=new XMLHttpRequest,h=V();if(e){if(y(a,"function"))m=b,b=a,a=null;else if(y(a,"object")){var d=[],f;for(f in a)a.hasOwnProperty(f)&&d.push(encodeURIComponent(f)+"="+encodeURIComponent(a[f]));a=d.join("&")}e.open(a?"POST":"GET",c,!0);a&&(e.setRequestHeader("X-Requested-With","XMLHttpRequest"),
e.setRequestHeader("Content-type","application/x-www-form-urlencoded"));b&&(k.once("snap.ajax."+h+".0",b),k.once("snap.ajax."+h+".200",b),k.once("snap.ajax."+h+".304",b));e.onreadystatechange=function(){4==e.readyState&&k("snap.ajax."+h+"."+e.status,m,e)};if(4==e.readyState)return e;e.send(a);return e}};a.load=function(c,b,m){a.ajax(c,function(c){c=a.parse(c.responseText);m?b.call(m,c):b(c)})};a.getElementByPoint=function(c,a){var b,m,e=G.doc.elementFromPoint(c,a);if(G.win.opera&&"svg"==e.tagName){b=
e;m=b.getBoundingClientRect();b=b.ownerDocument;var h=b.body,d=b.documentElement;b=m.top+(g.win.pageYOffset||d.scrollTop||h.scrollTop)-(d.clientTop||h.clientTop||0);m=m.left+(g.win.pageXOffset||d.scrollLeft||h.scrollLeft)-(d.clientLeft||h.clientLeft||0);h=e.createSVGRect();h.x=c-m;h.y=a-b;h.width=h.height=1;b=e.getIntersectionList(h,null);b.length&&(e=b[b.length-1])}return e?x(e):null};a.plugin=function(c){c(a,e,s,G,l)};return G.win.Snap=a}();C.plugin(function(a,k,y,M,A){function w(a,d,f,b,q,e){null==
d&&"[object SVGMatrix]"==z.call(a)?(this.a=a.a,this.b=a.b,this.c=a.c,this.d=a.d,this.e=a.e,this.f=a.f):null!=a?(this.a=+a,this.b=+d,this.c=+f,this.d=+b,this.e=+q,this.f=+e):(this.a=1,this.c=this.b=0,this.d=1,this.f=this.e=0)}var z=Object.prototype.toString,d=String,f=Math;(function(n){function k(a){return a[0]*a[0]+a[1]*a[1]}function p(a){var d=f.sqrt(k(a));a[0]&&(a[0]/=d);a[1]&&(a[1]/=d)}n.add=function(a,d,e,f,n,p){var k=[[],[],[] ],u=[[this.a,this.c,this.e],[this.b,this.d,this.f],[0,0,1] ];d=[[a,
e,n],[d,f,p],[0,0,1] ];a&&a instanceof w&&(d=[[a.a,a.c,a.e],[a.b,a.d,a.f],[0,0,1] ]);for(a=0;3>a;a++)for(e=0;3>e;e++){for(f=n=0;3>f;f++)n+=u[a][f]*d[f][e];k[a][e]=n}this.a=k[0][0];this.b=k[1][0];this.c=k[0][1];this.d=k[1][1];this.e=k[0][2];this.f=k[1][2];return this};n.invert=function(){var a=this.a*this.d-this.b*this.c;return new w(this.d/a,-this.b/a,-this.c/a,this.a/a,(this.c*this.f-this.d*this.e)/a,(this.b*this.e-this.a*this.f)/a)};n.clone=function(){return new w(this.a,this.b,this.c,this.d,this.e,
this.f)};n.translate=function(a,d){return this.add(1,0,0,1,a,d)};n.scale=function(a,d,e,f){null==d&&(d=a);(e||f)&&this.add(1,0,0,1,e,f);this.add(a,0,0,d,0,0);(e||f)&&this.add(1,0,0,1,-e,-f);return this};n.rotate=function(b,d,e){b=a.rad(b);d=d||0;e=e||0;var l=+f.cos(b).toFixed(9);b=+f.sin(b).toFixed(9);this.add(l,b,-b,l,d,e);return this.add(1,0,0,1,-d,-e)};n.x=function(a,d){return a*this.a+d*this.c+this.e};n.y=function(a,d){return a*this.b+d*this.d+this.f};n.get=function(a){return+this[d.fromCharCode(97+
a)].toFixed(4)};n.toString=function(){return"matrix("+[this.get(0),this.get(1),this.get(2),this.get(3),this.get(4),this.get(5)].join()+")"};n.offset=function(){return[this.e.toFixed(4),this.f.toFixed(4)]};n.determinant=function(){return this.a*this.d-this.b*this.c};n.split=function(){var b={};b.dx=this.e;b.dy=this.f;var d=[[this.a,this.c],[this.b,this.d] ];b.scalex=f.sqrt(k(d[0]));p(d[0]);b.shear=d[0][0]*d[1][0]+d[0][1]*d[1][1];d[1]=[d[1][0]-d[0][0]*b.shear,d[1][1]-d[0][1]*b.shear];b.scaley=f.sqrt(k(d[1]));
p(d[1]);b.shear/=b.scaley;0>this.determinant()&&(b.scalex=-b.scalex);var e=-d[0][1],d=d[1][1];0>d?(b.rotate=a.deg(f.acos(d)),0>e&&(b.rotate=360-b.rotate)):b.rotate=a.deg(f.asin(e));b.isSimple=!+b.shear.toFixed(9)&&(b.scalex.toFixed(9)==b.scaley.toFixed(9)||!b.rotate);b.isSuperSimple=!+b.shear.toFixed(9)&&b.scalex.toFixed(9)==b.scaley.toFixed(9)&&!b.rotate;b.noRotation=!+b.shear.toFixed(9)&&!b.rotate;return b};n.toTransformString=function(a){a=a||this.split();if(+a.shear.toFixed(9))return"m"+[this.get(0),
this.get(1),this.get(2),this.get(3),this.get(4),this.get(5)];a.scalex=+a.scalex.toFixed(4);a.scaley=+a.scaley.toFixed(4);a.rotate=+a.rotate.toFixed(4);return(a.dx||a.dy?"t"+[+a.dx.toFixed(4),+a.dy.toFixed(4)]:"")+(1!=a.scalex||1!=a.scaley?"s"+[a.scalex,a.scaley,0,0]:"")+(a.rotate?"r"+[+a.rotate.toFixed(4),0,0]:"")}})(w.prototype);a.Matrix=w;a.matrix=function(a,d,f,b,k,e){return new w(a,d,f,b,k,e)}});C.plugin(function(a,v,y,M,A){function w(h){return function(d){k.stop();d instanceof A&&1==d.node.childNodes.length&&
("radialGradient"==d.node.firstChild.tagName||"linearGradient"==d.node.firstChild.tagName||"pattern"==d.node.firstChild.tagName)&&(d=d.node.firstChild,b(this).appendChild(d),d=u(d));if(d instanceof v)if("radialGradient"==d.type||"linearGradient"==d.type||"pattern"==d.type){d.node.id||e(d.node,{id:d.id});var f=l(d.node.id)}else f=d.attr(h);else f=a.color(d),f.error?(f=a(b(this).ownerSVGElement).gradient(d))?(f.node.id||e(f.node,{id:f.id}),f=l(f.node.id)):f=d:f=r(f);d={};d[h]=f;e(this.node,d);this.node.style[h]=
x}}function z(a){k.stop();a==+a&&(a+="px");this.node.style.fontSize=a}function d(a){var b=[];a=a.childNodes;for(var e=0,f=a.length;e<f;e++){var l=a[e];3==l.nodeType&&b.push(l.nodeValue);"tspan"==l.tagName&&(1==l.childNodes.length&&3==l.firstChild.nodeType?b.push(l.firstChild.nodeValue):b.push(d(l)))}return b}function f(){k.stop();return this.node.style.fontSize}var n=a._.make,u=a._.wrap,p=a.is,b=a._.getSomeDefs,q=/^url\(#?([^)]+)\)$/,e=a._.$,l=a.url,r=String,s=a._.separator,x="";k.on("snap.util.attr.mask",
function(a){if(a instanceof v||a instanceof A){k.stop();a instanceof A&&1==a.node.childNodes.length&&(a=a.node.firstChild,b(this).appendChild(a),a=u(a));if("mask"==a.type)var d=a;else d=n("mask",b(this)),d.node.appendChild(a.node);!d.node.id&&e(d.node,{id:d.id});e(this.node,{mask:l(d.id)})}});(function(a){k.on("snap.util.attr.clip",a);k.on("snap.util.attr.clip-path",a);k.on("snap.util.attr.clipPath",a)})(function(a){if(a instanceof v||a instanceof A){k.stop();if("clipPath"==a.type)var d=a;else d=
n("clipPath",b(this)),d.node.appendChild(a.node),!d.node.id&&e(d.node,{id:d.id});e(this.node,{"clip-path":l(d.id)})}});k.on("snap.util.attr.fill",w("fill"));k.on("snap.util.attr.stroke",w("stroke"));var G=/^([lr])(?:\(([^)]*)\))?(.*)$/i;k.on("snap.util.grad.parse",function(a){a=r(a);var b=a.match(G);if(!b)return null;a=b[1];var e=b[2],b=b[3],e=e.split(/\s*,\s*/).map(function(a){return+a==a?+a:a});1==e.length&&0==e[0]&&(e=[]);b=b.split("-");b=b.map(function(a){a=a.split(":");var b={color:a[0]};a[1]&&
(b.offset=parseFloat(a[1]));return b});return{type:a,params:e,stops:b}});k.on("snap.util.attr.d",function(b){k.stop();p(b,"array")&&p(b[0],"array")&&(b=a.path.toString.call(b));b=r(b);b.match(/[ruo]/i)&&(b=a.path.toAbsolute(b));e(this.node,{d:b})})(-1);k.on("snap.util.attr.#text",function(a){k.stop();a=r(a);for(a=M.doc.createTextNode(a);this.node.firstChild;)this.node.removeChild(this.node.firstChild);this.node.appendChild(a)})(-1);k.on("snap.util.attr.path",function(a){k.stop();this.attr({d:a})})(-1);
k.on("snap.util.attr.class",function(a){k.stop();this.node.className.baseVal=a})(-1);k.on("snap.util.attr.viewBox",function(a){a=p(a,"object")&&"x"in a?[a.x,a.y,a.width,a.height].join(" "):p(a,"array")?a.join(" "):a;e(this.node,{viewBox:a});k.stop()})(-1);k.on("snap.util.attr.transform",function(a){this.transform(a);k.stop()})(-1);k.on("snap.util.attr.r",function(a){"rect"==this.type&&(k.stop(),e(this.node,{rx:a,ry:a}))})(-1);k.on("snap.util.attr.textpath",function(a){k.stop();if("text"==this.type){var d,
f;if(!a&&this.textPath){for(a=this.textPath;a.node.firstChild;)this.node.appendChild(a.node.firstChild);a.remove();delete this.textPath}else if(p(a,"string")?(d=b(this),a=u(d.parentNode).path(a),d.appendChild(a.node),d=a.id,a.attr({id:d})):(a=u(a),a instanceof v&&(d=a.attr("id"),d||(d=a.id,a.attr({id:d})))),d)if(a=this.textPath,f=this.node,a)a.attr({"xlink:href":"#"+d});else{for(a=e("textPath",{"xlink:href":"#"+d});f.firstChild;)a.appendChild(f.firstChild);f.appendChild(a);this.textPath=u(a)}}})(-1);
k.on("snap.util.attr.text",function(a){if("text"==this.type){for(var b=this.node,d=function(a){var b=e("tspan");if(p(a,"array"))for(var f=0;f<a.length;f++)b.appendChild(d(a[f]));else b.appendChild(M.doc.createTextNode(a));b.normalize&&b.normalize();return b};b.firstChild;)b.removeChild(b.firstChild);for(a=d(a);a.firstChild;)b.appendChild(a.firstChild)}k.stop()})(-1);k.on("snap.util.attr.fontSize",z)(-1);k.on("snap.util.attr.font-size",z)(-1);k.on("snap.util.getattr.transform",function(){k.stop();
return this.transform()})(-1);k.on("snap.util.getattr.textpath",function(){k.stop();return this.textPath})(-1);(function(){function b(d){return function(){k.stop();var b=M.doc.defaultView.getComputedStyle(this.node,null).getPropertyValue("marker-"+d);return"none"==b?b:a(M.doc.getElementById(b.match(q)[1]))}}function d(a){return function(b){k.stop();var d="marker"+a.charAt(0).toUpperCase()+a.substring(1);if(""==b||!b)this.node.style[d]="none";else if("marker"==b.type){var f=b.node.id;f||e(b.node,{id:b.id});
this.node.style[d]=l(f)}}}k.on("snap.util.getattr.marker-end",b("end"))(-1);k.on("snap.util.getattr.markerEnd",b("end"))(-1);k.on("snap.util.getattr.marker-start",b("start"))(-1);k.on("snap.util.getattr.markerStart",b("start"))(-1);k.on("snap.util.getattr.marker-mid",b("mid"))(-1);k.on("snap.util.getattr.markerMid",b("mid"))(-1);k.on("snap.util.attr.marker-end",d("end"))(-1);k.on("snap.util.attr.markerEnd",d("end"))(-1);k.on("snap.util.attr.marker-start",d("start"))(-1);k.on("snap.util.attr.markerStart",
d("start"))(-1);k.on("snap.util.attr.marker-mid",d("mid"))(-1);k.on("snap.util.attr.markerMid",d("mid"))(-1)})();k.on("snap.util.getattr.r",function(){if("rect"==this.type&&e(this.node,"rx")==e(this.node,"ry"))return k.stop(),e(this.node,"rx")})(-1);k.on("snap.util.getattr.text",function(){if("text"==this.type||"tspan"==this.type){k.stop();var a=d(this.node);return 1==a.length?a[0]:a}})(-1);k.on("snap.util.getattr.#text",function(){return this.node.textContent})(-1);k.on("snap.util.getattr.viewBox",
function(){k.stop();var b=e(this.node,"viewBox");if(b)return b=b.split(s),a._.box(+b[0],+b[1],+b[2],+b[3])})(-1);k.on("snap.util.getattr.points",function(){var a=e(this.node,"points");k.stop();if(a)return a.split(s)})(-1);k.on("snap.util.getattr.path",function(){var a=e(this.node,"d");k.stop();return a})(-1);k.on("snap.util.getattr.class",function(){return this.node.className.baseVal})(-1);k.on("snap.util.getattr.fontSize",f)(-1);k.on("snap.util.getattr.font-size",f)(-1)});C.plugin(function(a,v,y,
M,A){function w(a){return a}function z(a){return function(b){return+b.toFixed(3)+a}}var d={"+":function(a,b){return a+b},"-":function(a,b){return a-b},"/":function(a,b){return a/b},"*":function(a,b){return a*b}},f=String,n=/[a-z]+$/i,u=/^\s*([+\-\/*])\s*=\s*([\d.eE+\-]+)\s*([^\d\s]+)?\s*$/;k.on("snap.util.attr",function(a){if(a=f(a).match(u)){var b=k.nt(),b=b.substring(b.lastIndexOf(".")+1),q=this.attr(b),e={};k.stop();var l=a[3]||"",r=q.match(n),s=d[a[1] ];r&&r==l?a=s(parseFloat(q),+a[2]):(q=this.asPX(b),
a=s(this.asPX(b),this.asPX(b,a[2]+l)));isNaN(q)||isNaN(a)||(e[b]=a,this.attr(e))}})(-10);k.on("snap.util.equal",function(a,b){var q=f(this.attr(a)||""),e=f(b).match(u);if(e){k.stop();var l=e[3]||"",r=q.match(n),s=d[e[1] ];if(r&&r==l)return{from:parseFloat(q),to:s(parseFloat(q),+e[2]),f:z(r)};q=this.asPX(a);return{from:q,to:s(q,this.asPX(a,e[2]+l)),f:w}}})(-10)});C.plugin(function(a,v,y,M,A){var w=y.prototype,z=a.is;w.rect=function(a,d,k,p,b,q){var e;null==q&&(q=b);z(a,"object")&&"[object Object]"==
a?e=a:null!=a&&(e={x:a,y:d,width:k,height:p},null!=b&&(e.rx=b,e.ry=q));return this.el("rect",e)};w.circle=function(a,d,k){var p;z(a,"object")&&"[object Object]"==a?p=a:null!=a&&(p={cx:a,cy:d,r:k});return this.el("circle",p)};var d=function(){function a(){this.parentNode.removeChild(this)}return function(d,k){var p=M.doc.createElement("img"),b=M.doc.body;p.style.cssText="position:absolute;left:-9999em;top:-9999em";p.onload=function(){k.call(p);p.onload=p.onerror=null;b.removeChild(p)};p.onerror=a;
b.appendChild(p);p.src=d}}();w.image=function(f,n,k,p,b){var q=this.el("image");if(z(f,"object")&&"src"in f)q.attr(f);else if(null!=f){var e={"xlink:href":f,preserveAspectRatio:"none"};null!=n&&null!=k&&(e.x=n,e.y=k);null!=p&&null!=b?(e.width=p,e.height=b):d(f,function(){a._.$(q.node,{width:this.offsetWidth,height:this.offsetHeight})});a._.$(q.node,e)}return q};w.ellipse=function(a,d,k,p){var b;z(a,"object")&&"[object Object]"==a?b=a:null!=a&&(b={cx:a,cy:d,rx:k,ry:p});return this.el("ellipse",b)};
w.path=function(a){var d;z(a,"object")&&!z(a,"array")?d=a:a&&(d={d:a});return this.el("path",d)};w.group=w.g=function(a){var d=this.el("g");1==arguments.length&&a&&!a.type?d.attr(a):arguments.length&&d.add(Array.prototype.slice.call(arguments,0));return d};w.svg=function(a,d,k,p,b,q,e,l){var r={};z(a,"object")&&null==d?r=a:(null!=a&&(r.x=a),null!=d&&(r.y=d),null!=k&&(r.width=k),null!=p&&(r.height=p),null!=b&&null!=q&&null!=e&&null!=l&&(r.viewBox=[b,q,e,l]));return this.el("svg",r)};w.mask=function(a){var d=
this.el("mask");1==arguments.length&&a&&!a.type?d.attr(a):arguments.length&&d.add(Array.prototype.slice.call(arguments,0));return d};w.ptrn=function(a,d,k,p,b,q,e,l){if(z(a,"object"))var r=a;else arguments.length?(r={},null!=a&&(r.x=a),null!=d&&(r.y=d),null!=k&&(r.width=k),null!=p&&(r.height=p),null!=b&&null!=q&&null!=e&&null!=l&&(r.viewBox=[b,q,e,l])):r={patternUnits:"userSpaceOnUse"};return this.el("pattern",r)};w.use=function(a){return null!=a?(make("use",this.node),a instanceof v&&(a.attr("id")||
a.attr({id:ID()}),a=a.attr("id")),this.el("use",{"xlink:href":a})):v.prototype.use.call(this)};w.text=function(a,d,k){var p={};z(a,"object")?p=a:null!=a&&(p={x:a,y:d,text:k||""});return this.el("text",p)};w.line=function(a,d,k,p){var b={};z(a,"object")?b=a:null!=a&&(b={x1:a,x2:k,y1:d,y2:p});return this.el("line",b)};w.polyline=function(a){1<arguments.length&&(a=Array.prototype.slice.call(arguments,0));var d={};z(a,"object")&&!z(a,"array")?d=a:null!=a&&(d={points:a});return this.el("polyline",d)};
w.polygon=function(a){1<arguments.length&&(a=Array.prototype.slice.call(arguments,0));var d={};z(a,"object")&&!z(a,"array")?d=a:null!=a&&(d={points:a});return this.el("polygon",d)};(function(){function d(){return this.selectAll("stop")}function n(b,d){var f=e("stop"),k={offset:+d+"%"};b=a.color(b);k["stop-color"]=b.hex;1>b.opacity&&(k["stop-opacity"]=b.opacity);e(f,k);this.node.appendChild(f);return this}function u(){if("linearGradient"==this.type){var b=e(this.node,"x1")||0,d=e(this.node,"x2")||
1,f=e(this.node,"y1")||0,k=e(this.node,"y2")||0;return a._.box(b,f,math.abs(d-b),math.abs(k-f))}b=this.node.r||0;return a._.box((this.node.cx||0.5)-b,(this.node.cy||0.5)-b,2*b,2*b)}function p(a,d){function f(a,b){for(var d=(b-u)/(a-w),e=w;e<a;e++)h[e].offset=+(+u+d*(e-w)).toFixed(2);w=a;u=b}var n=k("snap.util.grad.parse",null,d).firstDefined(),p;if(!n)return null;n.params.unshift(a);p="l"==n.type.toLowerCase()?b.apply(0,n.params):q.apply(0,n.params);n.type!=n.type.toLowerCase()&&e(p.node,{gradientUnits:"userSpaceOnUse"});
var h=n.stops,n=h.length,u=0,w=0;n--;for(var v=0;v<n;v++)"offset"in h[v]&&f(v,h[v].offset);h[n].offset=h[n].offset||100;f(n,h[n].offset);for(v=0;v<=n;v++){var y=h[v];p.addStop(y.color,y.offset)}return p}function b(b,k,p,q,w){b=a._.make("linearGradient",b);b.stops=d;b.addStop=n;b.getBBox=u;null!=k&&e(b.node,{x1:k,y1:p,x2:q,y2:w});return b}function q(b,k,p,q,w,h){b=a._.make("radialGradient",b);b.stops=d;b.addStop=n;b.getBBox=u;null!=k&&e(b.node,{cx:k,cy:p,r:q});null!=w&&null!=h&&e(b.node,{fx:w,fy:h});
return b}var e=a._.$;w.gradient=function(a){return p(this.defs,a)};w.gradientLinear=function(a,d,e,f){return b(this.defs,a,d,e,f)};w.gradientRadial=function(a,b,d,e,f){return q(this.defs,a,b,d,e,f)};w.toString=function(){var b=this.node.ownerDocument,d=b.createDocumentFragment(),b=b.createElement("div"),e=this.node.cloneNode(!0);d.appendChild(b);b.appendChild(e);a._.$(e,{xmlns:"http://www.w3.org/2000/svg"});b=b.innerHTML;d.removeChild(d.firstChild);return b};w.clear=function(){for(var a=this.node.firstChild,
b;a;)b=a.nextSibling,"defs"!=a.tagName?a.parentNode.removeChild(a):w.clear.call({node:a}),a=b}})()});C.plugin(function(a,k,y,M){function A(a){var b=A.ps=A.ps||{};b[a]?b[a].sleep=100:b[a]={sleep:100};setTimeout(function(){for(var d in b)b[L](d)&&d!=a&&(b[d].sleep--,!b[d].sleep&&delete b[d])});return b[a]}function w(a,b,d,e){null==a&&(a=b=d=e=0);null==b&&(b=a.y,d=a.width,e=a.height,a=a.x);return{x:a,y:b,width:d,w:d,height:e,h:e,x2:a+d,y2:b+e,cx:a+d/2,cy:b+e/2,r1:F.min(d,e)/2,r2:F.max(d,e)/2,r0:F.sqrt(d*
d+e*e)/2,path:s(a,b,d,e),vb:[a,b,d,e].join(" ")}}function z(){return this.join(",").replace(N,"$1")}function d(a){a=C(a);a.toString=z;return a}function f(a,b,d,h,f,k,l,n,p){if(null==p)return e(a,b,d,h,f,k,l,n);if(0>p||e(a,b,d,h,f,k,l,n)<p)p=void 0;else{var q=0.5,O=1-q,s;for(s=e(a,b,d,h,f,k,l,n,O);0.01<Z(s-p);)q/=2,O+=(s<p?1:-1)*q,s=e(a,b,d,h,f,k,l,n,O);p=O}return u(a,b,d,h,f,k,l,n,p)}function n(b,d){function e(a){return+(+a).toFixed(3)}return a._.cacher(function(a,h,l){a instanceof k&&(a=a.attr("d"));
a=I(a);for(var n,p,D,q,O="",s={},c=0,t=0,r=a.length;t<r;t++){D=a[t];if("M"==D[0])n=+D[1],p=+D[2];else{q=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6]);if(c+q>h){if(d&&!s.start){n=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6],h-c);O+=["C"+e(n.start.x),e(n.start.y),e(n.m.x),e(n.m.y),e(n.x),e(n.y)];if(l)return O;s.start=O;O=["M"+e(n.x),e(n.y)+"C"+e(n.n.x),e(n.n.y),e(n.end.x),e(n.end.y),e(D[5]),e(D[6])].join();c+=q;n=+D[5];p=+D[6];continue}if(!b&&!d)return n=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6],h-c)}c+=q;n=+D[5];p=+D[6]}O+=
D.shift()+D}s.end=O;return n=b?c:d?s:u(n,p,D[0],D[1],D[2],D[3],D[4],D[5],1)},null,a._.clone)}function u(a,b,d,e,h,f,k,l,n){var p=1-n,q=ma(p,3),s=ma(p,2),c=n*n,t=c*n,r=q*a+3*s*n*d+3*p*n*n*h+t*k,q=q*b+3*s*n*e+3*p*n*n*f+t*l,s=a+2*n*(d-a)+c*(h-2*d+a),t=b+2*n*(e-b)+c*(f-2*e+b),x=d+2*n*(h-d)+c*(k-2*h+d),c=e+2*n*(f-e)+c*(l-2*f+e);a=p*a+n*d;b=p*b+n*e;h=p*h+n*k;f=p*f+n*l;l=90-180*F.atan2(s-x,t-c)/S;return{x:r,y:q,m:{x:s,y:t},n:{x:x,y:c},start:{x:a,y:b},end:{x:h,y:f},alpha:l}}function p(b,d,e,h,f,n,k,l){a.is(b,
"array")||(b=[b,d,e,h,f,n,k,l]);b=U.apply(null,b);return w(b.min.x,b.min.y,b.max.x-b.min.x,b.max.y-b.min.y)}function b(a,b,d){return b>=a.x&&b<=a.x+a.width&&d>=a.y&&d<=a.y+a.height}function q(a,d){a=w(a);d=w(d);return b(d,a.x,a.y)||b(d,a.x2,a.y)||b(d,a.x,a.y2)||b(d,a.x2,a.y2)||b(a,d.x,d.y)||b(a,d.x2,d.y)||b(a,d.x,d.y2)||b(a,d.x2,d.y2)||(a.x<d.x2&&a.x>d.x||d.x<a.x2&&d.x>a.x)&&(a.y<d.y2&&a.y>d.y||d.y<a.y2&&d.y>a.y)}function e(a,b,d,e,h,f,n,k,l){null==l&&(l=1);l=(1<l?1:0>l?0:l)/2;for(var p=[-0.1252,
0.1252,-0.3678,0.3678,-0.5873,0.5873,-0.7699,0.7699,-0.9041,0.9041,-0.9816,0.9816],q=[0.2491,0.2491,0.2335,0.2335,0.2032,0.2032,0.1601,0.1601,0.1069,0.1069,0.0472,0.0472],s=0,c=0;12>c;c++)var t=l*p[c]+l,r=t*(t*(-3*a+9*d-9*h+3*n)+6*a-12*d+6*h)-3*a+3*d,t=t*(t*(-3*b+9*e-9*f+3*k)+6*b-12*e+6*f)-3*b+3*e,s=s+q[c]*F.sqrt(r*r+t*t);return l*s}function l(a,b,d){a=I(a);b=I(b);for(var h,f,l,n,k,s,r,O,x,c,t=d?0:[],w=0,v=a.length;w<v;w++)if(x=a[w],"M"==x[0])h=k=x[1],f=s=x[2];else{"C"==x[0]?(x=[h,f].concat(x.slice(1)),
h=x[6],f=x[7]):(x=[h,f,h,f,k,s,k,s],h=k,f=s);for(var G=0,y=b.length;G<y;G++)if(c=b[G],"M"==c[0])l=r=c[1],n=O=c[2];else{"C"==c[0]?(c=[l,n].concat(c.slice(1)),l=c[6],n=c[7]):(c=[l,n,l,n,r,O,r,O],l=r,n=O);var z;var K=x,B=c;z=d;var H=p(K),J=p(B);if(q(H,J)){for(var H=e.apply(0,K),J=e.apply(0,B),H=~~(H/8),J=~~(J/8),U=[],A=[],F={},M=z?0:[],P=0;P<H+1;P++){var C=u.apply(0,K.concat(P/H));U.push({x:C.x,y:C.y,t:P/H})}for(P=0;P<J+1;P++)C=u.apply(0,B.concat(P/J)),A.push({x:C.x,y:C.y,t:P/J});for(P=0;P<H;P++)for(K=
0;K<J;K++){var Q=U[P],L=U[P+1],B=A[K],C=A[K+1],N=0.001>Z(L.x-Q.x)?"y":"x",S=0.001>Z(C.x-B.x)?"y":"x",R;R=Q.x;var Y=Q.y,V=L.x,ea=L.y,fa=B.x,ga=B.y,ha=C.x,ia=C.y;if(W(R,V)<X(fa,ha)||X(R,V)>W(fa,ha)||W(Y,ea)<X(ga,ia)||X(Y,ea)>W(ga,ia))R=void 0;else{var $=(R*ea-Y*V)*(fa-ha)-(R-V)*(fa*ia-ga*ha),aa=(R*ea-Y*V)*(ga-ia)-(Y-ea)*(fa*ia-ga*ha),ja=(R-V)*(ga-ia)-(Y-ea)*(fa-ha);if(ja){var $=$/ja,aa=aa/ja,ja=+$.toFixed(2),ba=+aa.toFixed(2);R=ja<+X(R,V).toFixed(2)||ja>+W(R,V).toFixed(2)||ja<+X(fa,ha).toFixed(2)||
ja>+W(fa,ha).toFixed(2)||ba<+X(Y,ea).toFixed(2)||ba>+W(Y,ea).toFixed(2)||ba<+X(ga,ia).toFixed(2)||ba>+W(ga,ia).toFixed(2)?void 0:{x:$,y:aa}}else R=void 0}R&&F[R.x.toFixed(4)]!=R.y.toFixed(4)&&(F[R.x.toFixed(4)]=R.y.toFixed(4),Q=Q.t+Z((R[N]-Q[N])/(L[N]-Q[N]))*(L.t-Q.t),B=B.t+Z((R[S]-B[S])/(C[S]-B[S]))*(C.t-B.t),0<=Q&&1>=Q&&0<=B&&1>=B&&(z?M++:M.push({x:R.x,y:R.y,t1:Q,t2:B})))}z=M}else z=z?0:[];if(d)t+=z;else{H=0;for(J=z.length;H<J;H++)z[H].segment1=w,z[H].segment2=G,z[H].bez1=x,z[H].bez2=c;t=t.concat(z)}}}return t}
function r(a){var b=A(a);if(b.bbox)return C(b.bbox);if(!a)return w();a=I(a);for(var d=0,e=0,h=[],f=[],l,n=0,k=a.length;n<k;n++)l=a[n],"M"==l[0]?(d=l[1],e=l[2],h.push(d),f.push(e)):(d=U(d,e,l[1],l[2],l[3],l[4],l[5],l[6]),h=h.concat(d.min.x,d.max.x),f=f.concat(d.min.y,d.max.y),d=l[5],e=l[6]);a=X.apply(0,h);l=X.apply(0,f);h=W.apply(0,h);f=W.apply(0,f);f=w(a,l,h-a,f-l);b.bbox=C(f);return f}function s(a,b,d,e,h){if(h)return[["M",+a+ +h,b],["l",d-2*h,0],["a",h,h,0,0,1,h,h],["l",0,e-2*h],["a",h,h,0,0,1,
-h,h],["l",2*h-d,0],["a",h,h,0,0,1,-h,-h],["l",0,2*h-e],["a",h,h,0,0,1,h,-h],["z"] ];a=[["M",a,b],["l",d,0],["l",0,e],["l",-d,0],["z"] ];a.toString=z;return a}function x(a,b,d,e,h){null==h&&null==e&&(e=d);a=+a;b=+b;d=+d;e=+e;if(null!=h){var f=Math.PI/180,l=a+d*Math.cos(-e*f);a+=d*Math.cos(-h*f);var n=b+d*Math.sin(-e*f);b+=d*Math.sin(-h*f);d=[["M",l,n],["A",d,d,0,+(180<h-e),0,a,b] ]}else d=[["M",a,b],["m",0,-e],["a",d,e,0,1,1,0,2*e],["a",d,e,0,1,1,0,-2*e],["z"] ];d.toString=z;return d}function G(b){var e=
A(b);if(e.abs)return d(e.abs);Q(b,"array")&&Q(b&&b[0],"array")||(b=a.parsePathString(b));if(!b||!b.length)return[["M",0,0] ];var h=[],f=0,l=0,n=0,k=0,p=0;"M"==b[0][0]&&(f=+b[0][1],l=+b[0][2],n=f,k=l,p++,h[0]=["M",f,l]);for(var q=3==b.length&&"M"==b[0][0]&&"R"==b[1][0].toUpperCase()&&"Z"==b[2][0].toUpperCase(),s,r,w=p,c=b.length;w<c;w++){h.push(s=[]);r=b[w];p=r[0];if(p!=p.toUpperCase())switch(s[0]=p.toUpperCase(),s[0]){case "A":s[1]=r[1];s[2]=r[2];s[3]=r[3];s[4]=r[4];s[5]=r[5];s[6]=+r[6]+f;s[7]=+r[7]+
l;break;case "V":s[1]=+r[1]+l;break;case "H":s[1]=+r[1]+f;break;case "R":for(var t=[f,l].concat(r.slice(1)),u=2,v=t.length;u<v;u++)t[u]=+t[u]+f,t[++u]=+t[u]+l;h.pop();h=h.concat(P(t,q));break;case "O":h.pop();t=x(f,l,r[1],r[2]);t.push(t[0]);h=h.concat(t);break;case "U":h.pop();h=h.concat(x(f,l,r[1],r[2],r[3]));s=["U"].concat(h[h.length-1].slice(-2));break;case "M":n=+r[1]+f,k=+r[2]+l;default:for(u=1,v=r.length;u<v;u++)s[u]=+r[u]+(u%2?f:l)}else if("R"==p)t=[f,l].concat(r.slice(1)),h.pop(),h=h.concat(P(t,
q)),s=["R"].concat(r.slice(-2));else if("O"==p)h.pop(),t=x(f,l,r[1],r[2]),t.push(t[0]),h=h.concat(t);else if("U"==p)h.pop(),h=h.concat(x(f,l,r[1],r[2],r[3])),s=["U"].concat(h[h.length-1].slice(-2));else for(t=0,u=r.length;t<u;t++)s[t]=r[t];p=p.toUpperCase();if("O"!=p)switch(s[0]){case "Z":f=+n;l=+k;break;case "H":f=s[1];break;case "V":l=s[1];break;case "M":n=s[s.length-2],k=s[s.length-1];default:f=s[s.length-2],l=s[s.length-1]}}h.toString=z;e.abs=d(h);return h}function h(a,b,d,e){return[a,b,d,e,d,
e]}function J(a,b,d,e,h,f){var l=1/3,n=2/3;return[l*a+n*d,l*b+n*e,l*h+n*d,l*f+n*e,h,f]}function K(b,d,e,h,f,l,n,k,p,s){var r=120*S/180,q=S/180*(+f||0),c=[],t,x=a._.cacher(function(a,b,c){var d=a*F.cos(c)-b*F.sin(c);a=a*F.sin(c)+b*F.cos(c);return{x:d,y:a}});if(s)v=s[0],t=s[1],l=s[2],u=s[3];else{t=x(b,d,-q);b=t.x;d=t.y;t=x(k,p,-q);k=t.x;p=t.y;F.cos(S/180*f);F.sin(S/180*f);t=(b-k)/2;v=(d-p)/2;u=t*t/(e*e)+v*v/(h*h);1<u&&(u=F.sqrt(u),e*=u,h*=u);var u=e*e,w=h*h,u=(l==n?-1:1)*F.sqrt(Z((u*w-u*v*v-w*t*t)/
(u*v*v+w*t*t)));l=u*e*v/h+(b+k)/2;var u=u*-h*t/e+(d+p)/2,v=F.asin(((d-u)/h).toFixed(9));t=F.asin(((p-u)/h).toFixed(9));v=b<l?S-v:v;t=k<l?S-t:t;0>v&&(v=2*S+v);0>t&&(t=2*S+t);n&&v>t&&(v-=2*S);!n&&t>v&&(t-=2*S)}if(Z(t-v)>r){var c=t,w=k,G=p;t=v+r*(n&&t>v?1:-1);k=l+e*F.cos(t);p=u+h*F.sin(t);c=K(k,p,e,h,f,0,n,w,G,[t,c,l,u])}l=t-v;f=F.cos(v);r=F.sin(v);n=F.cos(t);t=F.sin(t);l=F.tan(l/4);e=4/3*e*l;l*=4/3*h;h=[b,d];b=[b+e*r,d-l*f];d=[k+e*t,p-l*n];k=[k,p];b[0]=2*h[0]-b[0];b[1]=2*h[1]-b[1];if(s)return[b,d,k].concat(c);
c=[b,d,k].concat(c).join().split(",");s=[];k=0;for(p=c.length;k<p;k++)s[k]=k%2?x(c[k-1],c[k],q).y:x(c[k],c[k+1],q).x;return s}function U(a,b,d,e,h,f,l,k){for(var n=[],p=[[],[] ],s,r,c,t,q=0;2>q;++q)0==q?(r=6*a-12*d+6*h,s=-3*a+9*d-9*h+3*l,c=3*d-3*a):(r=6*b-12*e+6*f,s=-3*b+9*e-9*f+3*k,c=3*e-3*b),1E-12>Z(s)?1E-12>Z(r)||(s=-c/r,0<s&&1>s&&n.push(s)):(t=r*r-4*c*s,c=F.sqrt(t),0>t||(t=(-r+c)/(2*s),0<t&&1>t&&n.push(t),s=(-r-c)/(2*s),0<s&&1>s&&n.push(s)));for(r=q=n.length;q--;)s=n[q],c=1-s,p[0][q]=c*c*c*a+3*
c*c*s*d+3*c*s*s*h+s*s*s*l,p[1][q]=c*c*c*b+3*c*c*s*e+3*c*s*s*f+s*s*s*k;p[0][r]=a;p[1][r]=b;p[0][r+1]=l;p[1][r+1]=k;p[0].length=p[1].length=r+2;return{min:{x:X.apply(0,p[0]),y:X.apply(0,p[1])},max:{x:W.apply(0,p[0]),y:W.apply(0,p[1])}}}function I(a,b){var e=!b&&A(a);if(!b&&e.curve)return d(e.curve);var f=G(a),l=b&&G(b),n={x:0,y:0,bx:0,by:0,X:0,Y:0,qx:null,qy:null},k={x:0,y:0,bx:0,by:0,X:0,Y:0,qx:null,qy:null},p=function(a,b,c){if(!a)return["C",b.x,b.y,b.x,b.y,b.x,b.y];a[0]in{T:1,Q:1}||(b.qx=b.qy=null);
switch(a[0]){case "M":b.X=a[1];b.Y=a[2];break;case "A":a=["C"].concat(K.apply(0,[b.x,b.y].concat(a.slice(1))));break;case "S":"C"==c||"S"==c?(c=2*b.x-b.bx,b=2*b.y-b.by):(c=b.x,b=b.y);a=["C",c,b].concat(a.slice(1));break;case "T":"Q"==c||"T"==c?(b.qx=2*b.x-b.qx,b.qy=2*b.y-b.qy):(b.qx=b.x,b.qy=b.y);a=["C"].concat(J(b.x,b.y,b.qx,b.qy,a[1],a[2]));break;case "Q":b.qx=a[1];b.qy=a[2];a=["C"].concat(J(b.x,b.y,a[1],a[2],a[3],a[4]));break;case "L":a=["C"].concat(h(b.x,b.y,a[1],a[2]));break;case "H":a=["C"].concat(h(b.x,
b.y,a[1],b.y));break;case "V":a=["C"].concat(h(b.x,b.y,b.x,a[1]));break;case "Z":a=["C"].concat(h(b.x,b.y,b.X,b.Y))}return a},s=function(a,b){if(7<a[b].length){a[b].shift();for(var c=a[b];c.length;)q[b]="A",l&&(u[b]="A"),a.splice(b++,0,["C"].concat(c.splice(0,6)));a.splice(b,1);v=W(f.length,l&&l.length||0)}},r=function(a,b,c,d,e){a&&b&&"M"==a[e][0]&&"M"!=b[e][0]&&(b.splice(e,0,["M",d.x,d.y]),c.bx=0,c.by=0,c.x=a[e][1],c.y=a[e][2],v=W(f.length,l&&l.length||0))},q=[],u=[],c="",t="",x=0,v=W(f.length,
l&&l.length||0);for(;x<v;x++){f[x]&&(c=f[x][0]);"C"!=c&&(q[x]=c,x&&(t=q[x-1]));f[x]=p(f[x],n,t);"A"!=q[x]&&"C"==c&&(q[x]="C");s(f,x);l&&(l[x]&&(c=l[x][0]),"C"!=c&&(u[x]=c,x&&(t=u[x-1])),l[x]=p(l[x],k,t),"A"!=u[x]&&"C"==c&&(u[x]="C"),s(l,x));r(f,l,n,k,x);r(l,f,k,n,x);var w=f[x],z=l&&l[x],y=w.length,U=l&&z.length;n.x=w[y-2];n.y=w[y-1];n.bx=$(w[y-4])||n.x;n.by=$(w[y-3])||n.y;k.bx=l&&($(z[U-4])||k.x);k.by=l&&($(z[U-3])||k.y);k.x=l&&z[U-2];k.y=l&&z[U-1]}l||(e.curve=d(f));return l?[f,l]:f}function P(a,
b){for(var d=[],e=0,h=a.length;h-2*!b>e;e+=2){var f=[{x:+a[e-2],y:+a[e-1]},{x:+a[e],y:+a[e+1]},{x:+a[e+2],y:+a[e+3]},{x:+a[e+4],y:+a[e+5]}];b?e?h-4==e?f[3]={x:+a[0],y:+a[1]}:h-2==e&&(f[2]={x:+a[0],y:+a[1]},f[3]={x:+a[2],y:+a[3]}):f[0]={x:+a[h-2],y:+a[h-1]}:h-4==e?f[3]=f[2]:e||(f[0]={x:+a[e],y:+a[e+1]});d.push(["C",(-f[0].x+6*f[1].x+f[2].x)/6,(-f[0].y+6*f[1].y+f[2].y)/6,(f[1].x+6*f[2].x-f[3].x)/6,(f[1].y+6*f[2].y-f[3].y)/6,f[2].x,f[2].y])}return d}y=k.prototype;var Q=a.is,C=a._.clone,L="hasOwnProperty",
N=/,?([a-z]),?/gi,$=parseFloat,F=Math,S=F.PI,X=F.min,W=F.max,ma=F.pow,Z=F.abs;M=n(1);var na=n(),ba=n(0,1),V=a._unit2px;a.path=A;a.path.getTotalLength=M;a.path.getPointAtLength=na;a.path.getSubpath=function(a,b,d){if(1E-6>this.getTotalLength(a)-d)return ba(a,b).end;a=ba(a,d,1);return b?ba(a,b).end:a};y.getTotalLength=function(){if(this.node.getTotalLength)return this.node.getTotalLength()};y.getPointAtLength=function(a){return na(this.attr("d"),a)};y.getSubpath=function(b,d){return a.path.getSubpath(this.attr("d"),
b,d)};a._.box=w;a.path.findDotsAtSegment=u;a.path.bezierBBox=p;a.path.isPointInsideBBox=b;a.path.isBBoxIntersect=q;a.path.intersection=function(a,b){return l(a,b)};a.path.intersectionNumber=function(a,b){return l(a,b,1)};a.path.isPointInside=function(a,d,e){var h=r(a);return b(h,d,e)&&1==l(a,[["M",d,e],["H",h.x2+10] ],1)%2};a.path.getBBox=r;a.path.get={path:function(a){return a.attr("path")},circle:function(a){a=V(a);return x(a.cx,a.cy,a.r)},ellipse:function(a){a=V(a);return x(a.cx||0,a.cy||0,a.rx,
a.ry)},rect:function(a){a=V(a);return s(a.x||0,a.y||0,a.width,a.height,a.rx,a.ry)},image:function(a){a=V(a);return s(a.x||0,a.y||0,a.width,a.height)},line:function(a){return"M"+[a.attr("x1")||0,a.attr("y1")||0,a.attr("x2"),a.attr("y2")]},polyline:function(a){return"M"+a.attr("points")},polygon:function(a){return"M"+a.attr("points")+"z"},deflt:function(a){a=a.node.getBBox();return s(a.x,a.y,a.width,a.height)}};a.path.toRelative=function(b){var e=A(b),h=String.prototype.toLowerCase;if(e.rel)return d(e.rel);
a.is(b,"array")&&a.is(b&&b[0],"array")||(b=a.parsePathString(b));var f=[],l=0,n=0,k=0,p=0,s=0;"M"==b[0][0]&&(l=b[0][1],n=b[0][2],k=l,p=n,s++,f.push(["M",l,n]));for(var r=b.length;s<r;s++){var q=f[s]=[],x=b[s];if(x[0]!=h.call(x[0]))switch(q[0]=h.call(x[0]),q[0]){case "a":q[1]=x[1];q[2]=x[2];q[3]=x[3];q[4]=x[4];q[5]=x[5];q[6]=+(x[6]-l).toFixed(3);q[7]=+(x[7]-n).toFixed(3);break;case "v":q[1]=+(x[1]-n).toFixed(3);break;case "m":k=x[1],p=x[2];default:for(var c=1,t=x.length;c<t;c++)q[c]=+(x[c]-(c%2?l:
n)).toFixed(3)}else for(f[s]=[],"m"==x[0]&&(k=x[1]+l,p=x[2]+n),q=0,c=x.length;q<c;q++)f[s][q]=x[q];x=f[s].length;switch(f[s][0]){case "z":l=k;n=p;break;case "h":l+=+f[s][x-1];break;case "v":n+=+f[s][x-1];break;default:l+=+f[s][x-2],n+=+f[s][x-1]}}f.toString=z;e.rel=d(f);return f};a.path.toAbsolute=G;a.path.toCubic=I;a.path.map=function(a,b){if(!b)return a;var d,e,h,f,l,n,k;a=I(a);h=0;for(l=a.length;h<l;h++)for(k=a[h],f=1,n=k.length;f<n;f+=2)d=b.x(k[f],k[f+1]),e=b.y(k[f],k[f+1]),k[f]=d,k[f+1]=e;return a};
a.path.toString=z;a.path.clone=d});C.plugin(function(a,v,y,C){var A=Math.max,w=Math.min,z=function(a){this.items=[];this.bindings={};this.length=0;this.type="set";if(a)for(var f=0,n=a.length;f<n;f++)a[f]&&(this[this.items.length]=this.items[this.items.length]=a[f],this.length++)};v=z.prototype;v.push=function(){for(var a,f,n=0,k=arguments.length;n<k;n++)if(a=arguments[n])f=this.items.length,this[f]=this.items[f]=a,this.length++;return this};v.pop=function(){this.length&&delete this[this.length--];
return this.items.pop()};v.forEach=function(a,f){for(var n=0,k=this.items.length;n<k&&!1!==a.call(f,this.items[n],n);n++);return this};v.animate=function(d,f,n,u){"function"!=typeof n||n.length||(u=n,n=L.linear);d instanceof a._.Animation&&(u=d.callback,n=d.easing,f=n.dur,d=d.attr);var p=arguments;if(a.is(d,"array")&&a.is(p[p.length-1],"array"))var b=!0;var q,e=function(){q?this.b=q:q=this.b},l=0,r=u&&function(){l++==this.length&&u.call(this)};return this.forEach(function(a,l){k.once("snap.animcreated."+
a.id,e);b?p[l]&&a.animate.apply(a,p[l]):a.animate(d,f,n,r)})};v.remove=function(){for(;this.length;)this.pop().remove();return this};v.bind=function(a,f,k){var u={};if("function"==typeof f)this.bindings[a]=f;else{var p=k||a;this.bindings[a]=function(a){u[p]=a;f.attr(u)}}return this};v.attr=function(a){var f={},k;for(k in a)if(this.bindings[k])this.bindings[k](a[k]);else f[k]=a[k];a=0;for(k=this.items.length;a<k;a++)this.items[a].attr(f);return this};v.clear=function(){for(;this.length;)this.pop()};
v.splice=function(a,f,k){a=0>a?A(this.length+a,0):a;f=A(0,w(this.length-a,f));var u=[],p=[],b=[],q;for(q=2;q<arguments.length;q++)b.push(arguments[q]);for(q=0;q<f;q++)p.push(this[a+q]);for(;q<this.length-a;q++)u.push(this[a+q]);var e=b.length;for(q=0;q<e+u.length;q++)this.items[a+q]=this[a+q]=q<e?b[q]:u[q-e];for(q=this.items.length=this.length-=f-e;this[q];)delete this[q++];return new z(p)};v.exclude=function(a){for(var f=0,k=this.length;f<k;f++)if(this[f]==a)return this.splice(f,1),!0;return!1};
v.insertAfter=function(a){for(var f=this.items.length;f--;)this.items[f].insertAfter(a);return this};v.getBBox=function(){for(var a=[],f=[],k=[],u=[],p=this.items.length;p--;)if(!this.items[p].removed){var b=this.items[p].getBBox();a.push(b.x);f.push(b.y);k.push(b.x+b.width);u.push(b.y+b.height)}a=w.apply(0,a);f=w.apply(0,f);k=A.apply(0,k);u=A.apply(0,u);return{x:a,y:f,x2:k,y2:u,width:k-a,height:u-f,cx:a+(k-a)/2,cy:f+(u-f)/2}};v.clone=function(a){a=new z;for(var f=0,k=this.items.length;f<k;f++)a.push(this.items[f].clone());
return a};v.toString=function(){return"Snap\u2018s set"};v.type="set";a.set=function(){var a=new z;arguments.length&&a.push.apply(a,Array.prototype.slice.call(arguments,0));return a}});C.plugin(function(a,v,y,C){function A(a){var b=a[0];switch(b.toLowerCase()){case "t":return[b,0,0];case "m":return[b,1,0,0,1,0,0];case "r":return 4==a.length?[b,0,a[2],a[3] ]:[b,0];case "s":return 5==a.length?[b,1,1,a[3],a[4] ]:3==a.length?[b,1,1]:[b,1]}}function w(b,d,f){d=q(d).replace(/\.{3}|\u2026/g,b);b=a.parseTransformString(b)||
[];d=a.parseTransformString(d)||[];for(var k=Math.max(b.length,d.length),p=[],v=[],h=0,w,z,y,I;h<k;h++){y=b[h]||A(d[h]);I=d[h]||A(y);if(y[0]!=I[0]||"r"==y[0].toLowerCase()&&(y[2]!=I[2]||y[3]!=I[3])||"s"==y[0].toLowerCase()&&(y[3]!=I[3]||y[4]!=I[4])){b=a._.transform2matrix(b,f());d=a._.transform2matrix(d,f());p=[["m",b.a,b.b,b.c,b.d,b.e,b.f] ];v=[["m",d.a,d.b,d.c,d.d,d.e,d.f] ];break}p[h]=[];v[h]=[];w=0;for(z=Math.max(y.length,I.length);w<z;w++)w in y&&(p[h][w]=y[w]),w in I&&(v[h][w]=I[w])}return{from:u(p),
to:u(v),f:n(p)}}function z(a){return a}function d(a){return function(b){return+b.toFixed(3)+a}}function f(b){return a.rgb(b[0],b[1],b[2])}function n(a){var b=0,d,f,k,n,h,p,q=[];d=0;for(f=a.length;d<f;d++){h="[";p=['"'+a[d][0]+'"'];k=1;for(n=a[d].length;k<n;k++)p[k]="val["+b++ +"]";h+=p+"]";q[d]=h}return Function("val","return Snap.path.toString.call(["+q+"])")}function u(a){for(var b=[],d=0,f=a.length;d<f;d++)for(var k=1,n=a[d].length;k<n;k++)b.push(a[d][k]);return b}var p={},b=/[a-z]+$/i,q=String;
p.stroke=p.fill="colour";v.prototype.equal=function(a,b){return k("snap.util.equal",this,a,b).firstDefined()};k.on("snap.util.equal",function(e,k){var r,s;r=q(this.attr(e)||"");var x=this;if(r==+r&&k==+k)return{from:+r,to:+k,f:z};if("colour"==p[e])return r=a.color(r),s=a.color(k),{from:[r.r,r.g,r.b,r.opacity],to:[s.r,s.g,s.b,s.opacity],f:f};if("transform"==e||"gradientTransform"==e||"patternTransform"==e)return k instanceof a.Matrix&&(k=k.toTransformString()),a._.rgTransform.test(k)||(k=a._.svgTransform2string(k)),
w(r,k,function(){return x.getBBox(1)});if("d"==e||"path"==e)return r=a.path.toCubic(r,k),{from:u(r[0]),to:u(r[1]),f:n(r[0])};if("points"==e)return r=q(r).split(a._.separator),s=q(k).split(a._.separator),{from:r,to:s,f:function(a){return a}};aUnit=r.match(b);s=q(k).match(b);return aUnit&&aUnit==s?{from:parseFloat(r),to:parseFloat(k),f:d(aUnit)}:{from:this.asPX(e),to:this.asPX(e,k),f:z}})});C.plugin(function(a,v,y,C){var A=v.prototype,w="createTouch"in C.doc;v="click dblclick mousedown mousemove mouseout mouseover mouseup touchstart touchmove touchend touchcancel".split(" ");
var z={mousedown:"touchstart",mousemove:"touchmove",mouseup:"touchend"},d=function(a,b){var d="y"==a?"scrollTop":"scrollLeft",e=b&&b.node?b.node.ownerDocument:C.doc;return e[d in e.documentElement?"documentElement":"body"][d]},f=function(){this.returnValue=!1},n=function(){return this.originalEvent.preventDefault()},u=function(){this.cancelBubble=!0},p=function(){return this.originalEvent.stopPropagation()},b=function(){if(C.doc.addEventListener)return function(a,b,e,f){var k=w&&z[b]?z[b]:b,l=function(k){var l=
d("y",f),q=d("x",f);if(w&&z.hasOwnProperty(b))for(var r=0,u=k.targetTouches&&k.targetTouches.length;r<u;r++)if(k.targetTouches[r].target==a||a.contains(k.targetTouches[r].target)){u=k;k=k.targetTouches[r];k.originalEvent=u;k.preventDefault=n;k.stopPropagation=p;break}return e.call(f,k,k.clientX+q,k.clientY+l)};b!==k&&a.addEventListener(b,l,!1);a.addEventListener(k,l,!1);return function(){b!==k&&a.removeEventListener(b,l,!1);a.removeEventListener(k,l,!1);return!0}};if(C.doc.attachEvent)return function(a,
b,e,h){var k=function(a){a=a||h.node.ownerDocument.window.event;var b=d("y",h),k=d("x",h),k=a.clientX+k,b=a.clientY+b;a.preventDefault=a.preventDefault||f;a.stopPropagation=a.stopPropagation||u;return e.call(h,a,k,b)};a.attachEvent("on"+b,k);return function(){a.detachEvent("on"+b,k);return!0}}}(),q=[],e=function(a){for(var b=a.clientX,e=a.clientY,f=d("y"),l=d("x"),n,p=q.length;p--;){n=q[p];if(w)for(var r=a.touches&&a.touches.length,u;r--;){if(u=a.touches[r],u.identifier==n.el._drag.id||n.el.node.contains(u.target)){b=
u.clientX;e=u.clientY;(a.originalEvent?a.originalEvent:a).preventDefault();break}}else a.preventDefault();b+=l;e+=f;k("snap.drag.move."+n.el.id,n.move_scope||n.el,b-n.el._drag.x,e-n.el._drag.y,b,e,a)}},l=function(b){a.unmousemove(e).unmouseup(l);for(var d=q.length,f;d--;)f=q[d],f.el._drag={},k("snap.drag.end."+f.el.id,f.end_scope||f.start_scope||f.move_scope||f.el,b);q=[]};for(y=v.length;y--;)(function(d){a[d]=A[d]=function(e,f){a.is(e,"function")&&(this.events=this.events||[],this.events.push({name:d,
f:e,unbind:b(this.node||document,d,e,f||this)}));return this};a["un"+d]=A["un"+d]=function(a){for(var b=this.events||[],e=b.length;e--;)if(b[e].name==d&&(b[e].f==a||!a)){b[e].unbind();b.splice(e,1);!b.length&&delete this.events;break}return this}})(v[y]);A.hover=function(a,b,d,e){return this.mouseover(a,d).mouseout(b,e||d)};A.unhover=function(a,b){return this.unmouseover(a).unmouseout(b)};var r=[];A.drag=function(b,d,f,h,n,p){function u(r,v,w){(r.originalEvent||r).preventDefault();this._drag.x=v;
this._drag.y=w;this._drag.id=r.identifier;!q.length&&a.mousemove(e).mouseup(l);q.push({el:this,move_scope:h,start_scope:n,end_scope:p});d&&k.on("snap.drag.start."+this.id,d);b&&k.on("snap.drag.move."+this.id,b);f&&k.on("snap.drag.end."+this.id,f);k("snap.drag.start."+this.id,n||h||this,v,w,r)}if(!arguments.length){var v;return this.drag(function(a,b){this.attr({transform:v+(v?"T":"t")+[a,b]})},function(){v=this.transform().local})}this._drag={};r.push({el:this,start:u});this.mousedown(u);return this};
A.undrag=function(){for(var b=r.length;b--;)r[b].el==this&&(this.unmousedown(r[b].start),r.splice(b,1),k.unbind("snap.drag.*."+this.id));!r.length&&a.unmousemove(e).unmouseup(l);return this}});C.plugin(function(a,v,y,C){y=y.prototype;var A=/^\s*url\((.+)\)/,w=String,z=a._.$;a.filter={};y.filter=function(d){var f=this;"svg"!=f.type&&(f=f.paper);d=a.parse(w(d));var k=a._.id(),u=z("filter");z(u,{id:k,filterUnits:"userSpaceOnUse"});u.appendChild(d.node);f.defs.appendChild(u);return new v(u)};k.on("snap.util.getattr.filter",
function(){k.stop();var d=z(this.node,"filter");if(d)return(d=w(d).match(A))&&a.select(d[1])});k.on("snap.util.attr.filter",function(d){if(d instanceof v&&"filter"==d.type){k.stop();var f=d.node.id;f||(z(d.node,{id:d.id}),f=d.id);z(this.node,{filter:a.url(f)})}d&&"none"!=d||(k.stop(),this.node.removeAttribute("filter"))});a.filter.blur=function(d,f){null==d&&(d=2);return a.format('<feGaussianBlur stdDeviation="{def}"/>',{def:null==f?d:[d,f]})};a.filter.blur.toString=function(){return this()};a.filter.shadow=
function(d,f,k,u,p){"string"==typeof k&&(p=u=k,k=4);"string"!=typeof u&&(p=u,u="#000");null==k&&(k=4);null==p&&(p=1);null==d&&(d=0,f=2);null==f&&(f=d);u=a.color(u||"#000");return a.format('<feGaussianBlur in="SourceAlpha" stdDeviation="{blur}"/><feOffset dx="{dx}" dy="{dy}" result="offsetblur"/><feFlood flood-color="{color}"/><feComposite in2="offsetblur" operator="in"/><feComponentTransfer><feFuncA type="linear" slope="{opacity}"/></feComponentTransfer><feMerge><feMergeNode/><feMergeNode in="SourceGraphic"/></feMerge>',
{color:u,dx:d,dy:f,blur:k,opacity:p})};a.filter.shadow.toString=function(){return this()};a.filter.grayscale=function(d){null==d&&(d=1);return a.format('<feColorMatrix type="matrix" values="{a} {b} {c} 0 0 {d} {e} {f} 0 0 {g} {b} {h} 0 0 0 0 0 1 0"/>',{a:0.2126+0.7874*(1-d),b:0.7152-0.7152*(1-d),c:0.0722-0.0722*(1-d),d:0.2126-0.2126*(1-d),e:0.7152+0.2848*(1-d),f:0.0722-0.0722*(1-d),g:0.2126-0.2126*(1-d),h:0.0722+0.9278*(1-d)})};a.filter.grayscale.toString=function(){return this()};a.filter.sepia=
function(d){null==d&&(d=1);return a.format('<feColorMatrix type="matrix" values="{a} {b} {c} 0 0 {d} {e} {f} 0 0 {g} {h} {i} 0 0 0 0 0 1 0"/>',{a:0.393+0.607*(1-d),b:0.769-0.769*(1-d),c:0.189-0.189*(1-d),d:0.349-0.349*(1-d),e:0.686+0.314*(1-d),f:0.168-0.168*(1-d),g:0.272-0.272*(1-d),h:0.534-0.534*(1-d),i:0.131+0.869*(1-d)})};a.filter.sepia.toString=function(){return this()};a.filter.saturate=function(d){null==d&&(d=1);return a.format('<feColorMatrix type="saturate" values="{amount}"/>',{amount:1-
d})};a.filter.saturate.toString=function(){return this()};a.filter.hueRotate=function(d){return a.format('<feColorMatrix type="hueRotate" values="{angle}"/>',{angle:d||0})};a.filter.hueRotate.toString=function(){return this()};a.filter.invert=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="table" tableValues="{amount} {amount2}"/><feFuncG type="table" tableValues="{amount} {amount2}"/><feFuncB type="table" tableValues="{amount} {amount2}"/></feComponentTransfer>',{amount:d,
amount2:1-d})};a.filter.invert.toString=function(){return this()};a.filter.brightness=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="linear" slope="{amount}"/><feFuncG type="linear" slope="{amount}"/><feFuncB type="linear" slope="{amount}"/></feComponentTransfer>',{amount:d})};a.filter.brightness.toString=function(){return this()};a.filter.contrast=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="linear" slope="{amount}" intercept="{amount2}"/><feFuncG type="linear" slope="{amount}" intercept="{amount2}"/><feFuncB type="linear" slope="{amount}" intercept="{amount2}"/></feComponentTransfer>',
{amount:d,amount2:0.5-d/2})};a.filter.contrast.toString=function(){return this()}});return C});

]]> </script>
</svg>



<div class="markdown"><h3>Block properties</h3>
</div>

<pre class='language-julia'><code class='language-julia'>mat(X)</code></pre>
<pre id='var-hash993548' class='code-output documenter-example-output'>2×2 LuxurySparse.SDPermMatrix{ComplexF64, Int64, Vector{ComplexF64}, Vector{Int64}}:
 0.0+0.0im  1.0+0.0im
 1.0+0.0im  0.0+0.0im</pre>


<div class="markdown"><ol>
<li><p>Hermitian: <span class="tex">$A^\dagger &#61; A$</span></p>
</li>
</ol>
</div>

<pre class='language-julia'><code class='language-julia'>ishermitian(X)</code></pre>
<pre id='var-hash165882' class='code-output documenter-example-output'>true</pre>


<div class="markdown"><p>Only Hermitian operators are allowed in time evolution.</p>
</div>


<div class="markdown"><p class="tex">$$\texttt&#123;time\_evolve&#40;H, t&#41;&#125; &#61; e^&#123;-iHt&#125;$$</p>
<p>When the input operator is not Hermitian, the block constructor will throw you an error:</p>
</div>

<pre class='language-julia'><code class='language-julia'># time_evolve(ConstGate.S, 0.4)</code></pre>


<pre class='language-julia'><code class='language-julia'>mat(ConstGate.S)  # this is not hermitian</code></pre>
<pre id='var-hash103460' class='code-output documenter-example-output'>2×2 Diagonal{ComplexF64, Vector{ComplexF64}}:
 1.0+0.0im      ⋅    
     ⋅      0.0+1.0im</pre>


<div class="markdown"><p>Hermitian operators are not changed by conjugate transpose</p>
</div>

<pre class='language-julia'><code class='language-julia'>qft3'</code></pre>
<pre id='var-hash982563' class='code-output documenter-example-output'>nqubits: 3
chain
├─ chain
│  └─ put on (3)
│     └─ H
├─ chain
│  ├─ control(3)
│  │  └─ (2,) shift(-1.5707963267948966)
│  └─ put on (2)
│     └─ H
└─ chain
   ├─ control(3)
   │  └─ (1,) shift(-0.7853981633974483)
   ├─ control(2)
   │  └─ (1,) shift(-1.5707963267948966)
   └─ put on (1)
      └─ H
</pre>

<pre class='language-julia'><code class='language-julia'>vizcircuit(qft3')</code></pre>
<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:xlink="http://www.w3.org/1999/xlink"
     version="1.2"
     width="148.75mm" height="56.25mm" viewBox="0 0 148.75 56.25"
     stroke="none"
     fill="#000000"
     stroke-width="0.3"
     font-size="3.88"

     id="img-4069dc56">
<defs>
  <marker id="arrow" markerWidth="15" markerHeight="7" refX="5" refY="3.5" orient="auto" markerUnits="strokeWidth">
    <path d="M0,0 L15,3.5 L0,7 z" stroke="context-stroke" fill="context-stroke"/>
  </marker>
</defs>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-4069dc56-1">
  <g transform="translate(129.37,6.25)">
    <path d="M-5,-5 L5,-5 5,5 -5,5  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0" fill="#000000" id="img-4069dc56-2">
  <g transform="translate(108.12,25)">
    <circle cx="0" cy="0" r="0.83" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-4069dc56-3">
  <g transform="translate(108.12,6.25)">
    <path d="M-8.31,-5 L8.31,-5 8.31,5 -8.31,5  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0" fill="#000000" id="img-4069dc56-4">
  <g transform="translate(86.87,43.75)">
    <circle cx="0" cy="0" r="0.83" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-4069dc56-5">
  <g transform="translate(86.87,6.25)">
    <path d="M-8.31,-5 L8.31,-5 8.31,5 -8.31,5  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-4069dc56-6">
  <g transform="translate(65.62,25)">
    <path d="M-5,-5 L5,-5 5,5 -5,5  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0" fill="#000000" id="img-4069dc56-7">
  <g transform="translate(44.37,43.75)">
    <circle cx="0" cy="0" r="0.83" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-4069dc56-8">
  <g transform="translate(44.37,25)">
    <path d="M-8.31,-5 L8.31,-5 8.31,5 -8.31,5  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-4069dc56-9">
  <g transform="translate(23.12,43.75)">
    <path d="M-5,-5 L5,-5 5,5 -5,5  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-4069dc56-10">
  <g transform="translate(118.75,43.75)">
    <path fill="none" d="M21.25,0 L-21.25,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-4069dc56-11">
  <g transform="translate(129.37,25)">
    <path fill="none" d="M10.62,0 L-10.62,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-4069dc56-12">
  <g transform="translate(140,6.25)">
    <path fill="none" d="M0,0 L0,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-4069dc56-13">
  <g transform="translate(137.19,6.25)">
    <path fill="none" d="M-2.81,0 L2.81,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-4069dc56-14">
  <g transform="translate(121.56,6.25)">
    <path fill="none" d="M-2.81,0 L2.81,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-4069dc56-15">
  <g transform="translate(113.85,25)">
    <path fill="none" d="M-4.9,0 L4.9,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-4069dc56-16">
  <g transform="translate(102.4,25)">
    <path fill="none" d="M-4.9,0 L4.9,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-4069dc56-17">
  <g transform="translate(117.59,6.25)">
    <path fill="none" d="M-1.16,0 L1.16,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-4069dc56-18">
  <g transform="translate(98.66,6.25)">
    <path fill="none" d="M-1.16,0 L1.16,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-4069dc56-19">
  <g transform="translate(108.12,17.71)">
    <path fill="none" d="M0,6.46 L0,-6.46 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-4069dc56-20">
  <g transform="translate(86.87,25)">
    <path fill="none" d="M-10.62,0 L10.62,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-4069dc56-21">
  <g transform="translate(92.6,43.75)">
    <path fill="none" d="M-4.9,0 L4.9,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-4069dc56-22">
  <g transform="translate(70.52,43.75)">
    <path fill="none" d="M-15.52,0 L15.52,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-4069dc56-23">
  <g transform="translate(96.34,6.25)">
    <path fill="none" d="M-1.16,0 L1.16,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-4069dc56-24">
  <g transform="translate(45.53,6.25)">
    <path fill="none" d="M-33.03,0 L33.03,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-4069dc56-25">
  <g transform="translate(86.87,27.08)">
    <path fill="none" d="M0,15.83 L0,-15.83 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-4069dc56-26">
  <g transform="translate(73.44,25)">
    <path fill="none" d="M-2.81,0 L2.81,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-4069dc56-27">
  <g transform="translate(57.81,25)">
    <path fill="none" d="M-2.81,0 L2.81,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-4069dc56-28">
  <g transform="translate(50.1,43.75)">
    <path fill="none" d="M-4.9,0 L4.9,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-4069dc56-29">
  <g transform="translate(38.65,43.75)">
    <path fill="none" d="M-4.9,0 L4.9,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-4069dc56-30">
  <g transform="translate(53.84,25)">
    <path fill="none" d="M-1.16,0 L1.16,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-4069dc56-31">
  <g transform="translate(24.28,25)">
    <path fill="none" d="M-11.78,0 L11.78,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-4069dc56-32">
  <g transform="translate(44.37,36.46)">
    <path fill="none" d="M0,6.46 L0,-6.46 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-4069dc56-33">
  <g transform="translate(30.94,43.75)">
    <path fill="none" d="M-2.81,0 L2.81,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-4069dc56-34">
  <g transform="translate(15.31,43.75)">
    <path fill="none" d="M-2.81,0 L2.81,0 " class="primitive"/>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="5.64" id="img-4069dc56-35">
  <g transform="translate(129.37,6.25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">H</text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="5.64" id="img-4069dc56-36">
  <g transform="translate(108.12,25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em"></text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="3.53" id="img-4069dc56-37">
  <g transform="translate(108.12,6.25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">ϕ(-π/2)</text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="5.64" id="img-4069dc56-38">
  <g transform="translate(86.87,43.75)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em"></text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="3.53" id="img-4069dc56-39">
  <g transform="translate(86.87,6.25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">ϕ(-π/4)</text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="5.64" id="img-4069dc56-40">
  <g transform="translate(65.62,25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">H</text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="5.64" id="img-4069dc56-41">
  <g transform="translate(44.37,43.75)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em"></text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="3.53" id="img-4069dc56-42">
  <g transform="translate(44.37,25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">ϕ(-π/2)</text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="5.64" id="img-4069dc56-43">
  <g transform="translate(23.12,43.75)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">H</text>
    </g>
  </g>
</g>
<script> <![CDATA[
(function(N){var k=/[\.\/]/,L=/\s*,\s*/,C=function(a,d){return a-d},a,v,y={n:{}},M=function(){for(var a=0,d=this.length;a<d;a++)if("undefined"!=typeof this[a])return this[a]},A=function(){for(var a=this.length;--a;)if("undefined"!=typeof this[a])return this[a]},w=function(k,d){k=String(k);var f=v,n=Array.prototype.slice.call(arguments,2),u=w.listeners(k),p=0,b,q=[],e={},l=[],r=a;l.firstDefined=M;l.lastDefined=A;a=k;for(var s=v=0,x=u.length;s<x;s++)"zIndex"in u[s]&&(q.push(u[s].zIndex),0>u[s].zIndex&&
(e[u[s].zIndex]=u[s]));for(q.sort(C);0>q[p];)if(b=e[q[p++] ],l.push(b.apply(d,n)),v)return v=f,l;for(s=0;s<x;s++)if(b=u[s],"zIndex"in b)if(b.zIndex==q[p]){l.push(b.apply(d,n));if(v)break;do if(p++,(b=e[q[p] ])&&l.push(b.apply(d,n)),v)break;while(b)}else e[b.zIndex]=b;else if(l.push(b.apply(d,n)),v)break;v=f;a=r;return l};w._events=y;w.listeners=function(a){a=a.split(k);var d=y,f,n,u,p,b,q,e,l=[d],r=[];u=0;for(p=a.length;u<p;u++){e=[];b=0;for(q=l.length;b<q;b++)for(d=l[b].n,f=[d[a[u] ],d["*"] ],n=2;n--;)if(d=
f[n])e.push(d),r=r.concat(d.f||[]);l=e}return r};w.on=function(a,d){a=String(a);if("function"!=typeof d)return function(){};for(var f=a.split(L),n=0,u=f.length;n<u;n++)(function(a){a=a.split(k);for(var b=y,f,e=0,l=a.length;e<l;e++)b=b.n,b=b.hasOwnProperty(a[e])&&b[a[e] ]||(b[a[e] ]={n:{}});b.f=b.f||[];e=0;for(l=b.f.length;e<l;e++)if(b.f[e]==d){f=!0;break}!f&&b.f.push(d)})(f[n]);return function(a){+a==+a&&(d.zIndex=+a)}};w.f=function(a){var d=[].slice.call(arguments,1);return function(){w.apply(null,
[a,null].concat(d).concat([].slice.call(arguments,0)))}};w.stop=function(){v=1};w.nt=function(k){return k?(new RegExp("(?:\\.|\\/|^)"+k+"(?:\\.|\\/|$)")).test(a):a};w.nts=function(){return a.split(k)};w.off=w.unbind=function(a,d){if(a){var f=a.split(L);if(1<f.length)for(var n=0,u=f.length;n<u;n++)w.off(f[n],d);else{for(var f=a.split(k),p,b,q,e,l=[y],n=0,u=f.length;n<u;n++)for(e=0;e<l.length;e+=q.length-2){q=[e,1];p=l[e].n;if("*"!=f[n])p[f[n] ]&&q.push(p[f[n] ]);else for(b in p)p.hasOwnProperty(b)&&
q.push(p[b]);l.splice.apply(l,q)}n=0;for(u=l.length;n<u;n++)for(p=l[n];p.n;){if(d){if(p.f){e=0;for(f=p.f.length;e<f;e++)if(p.f[e]==d){p.f.splice(e,1);break}!p.f.length&&delete p.f}for(b in p.n)if(p.n.hasOwnProperty(b)&&p.n[b].f){q=p.n[b].f;e=0;for(f=q.length;e<f;e++)if(q[e]==d){q.splice(e,1);break}!q.length&&delete p.n[b].f}}else for(b in delete p.f,p.n)p.n.hasOwnProperty(b)&&p.n[b].f&&delete p.n[b].f;p=p.n}}}else w._events=y={n:{}}};w.once=function(a,d){var f=function(){w.unbind(a,f);return d.apply(this,
arguments)};return w.on(a,f)};w.version="0.4.2";w.toString=function(){return"You are running Eve 0.4.2"};"undefined"!=typeof module&&module.exports?module.exports=w:"function"===typeof define&&define.amd?define("eve",[],function(){return w}):N.eve=w})(this);
(function(N,k){"function"===typeof define&&define.amd?define("Snap.svg",["eve"],function(L){return k(N,L)}):k(N,N.eve)})(this,function(N,k){var L=function(a){var k={},y=N.requestAnimationFrame||N.webkitRequestAnimationFrame||N.mozRequestAnimationFrame||N.oRequestAnimationFrame||N.msRequestAnimationFrame||function(a){setTimeout(a,16)},M=Array.isArray||function(a){return a instanceof Array||"[object Array]"==Object.prototype.toString.call(a)},A=0,w="M"+(+new Date).toString(36),z=function(a){if(null==
a)return this.s;var b=this.s-a;this.b+=this.dur*b;this.B+=this.dur*b;this.s=a},d=function(a){if(null==a)return this.spd;this.spd=a},f=function(a){if(null==a)return this.dur;this.s=this.s*a/this.dur;this.dur=a},n=function(){delete k[this.id];this.update();a("mina.stop."+this.id,this)},u=function(){this.pdif||(delete k[this.id],this.update(),this.pdif=this.get()-this.b)},p=function(){this.pdif&&(this.b=this.get()-this.pdif,delete this.pdif,k[this.id]=this)},b=function(){var a;if(M(this.start)){a=[];
for(var b=0,e=this.start.length;b<e;b++)a[b]=+this.start[b]+(this.end[b]-this.start[b])*this.easing(this.s)}else a=+this.start+(this.end-this.start)*this.easing(this.s);this.set(a)},q=function(){var l=0,b;for(b in k)if(k.hasOwnProperty(b)){var e=k[b],f=e.get();l++;e.s=(f-e.b)/(e.dur/e.spd);1<=e.s&&(delete k[b],e.s=1,l--,function(b){setTimeout(function(){a("mina.finish."+b.id,b)})}(e));e.update()}l&&y(q)},e=function(a,r,s,x,G,h,J){a={id:w+(A++).toString(36),start:a,end:r,b:s,s:0,dur:x-s,spd:1,get:G,
set:h,easing:J||e.linear,status:z,speed:d,duration:f,stop:n,pause:u,resume:p,update:b};k[a.id]=a;r=0;for(var K in k)if(k.hasOwnProperty(K)&&(r++,2==r))break;1==r&&y(q);return a};e.time=Date.now||function(){return+new Date};e.getById=function(a){return k[a]||null};e.linear=function(a){return a};e.easeout=function(a){return Math.pow(a,1.7)};e.easein=function(a){return Math.pow(a,0.48)};e.easeinout=function(a){if(1==a)return 1;if(0==a)return 0;var b=0.48-a/1.04,e=Math.sqrt(0.1734+b*b);a=e-b;a=Math.pow(Math.abs(a),
1/3)*(0>a?-1:1);b=-e-b;b=Math.pow(Math.abs(b),1/3)*(0>b?-1:1);a=a+b+0.5;return 3*(1-a)*a*a+a*a*a};e.backin=function(a){return 1==a?1:a*a*(2.70158*a-1.70158)};e.backout=function(a){if(0==a)return 0;a-=1;return a*a*(2.70158*a+1.70158)+1};e.elastic=function(a){return a==!!a?a:Math.pow(2,-10*a)*Math.sin(2*(a-0.075)*Math.PI/0.3)+1};e.bounce=function(a){a<1/2.75?a*=7.5625*a:a<2/2.75?(a-=1.5/2.75,a=7.5625*a*a+0.75):a<2.5/2.75?(a-=2.25/2.75,a=7.5625*a*a+0.9375):(a-=2.625/2.75,a=7.5625*a*a+0.984375);return a};
return N.mina=e}("undefined"==typeof k?function(){}:k),C=function(){function a(c,t){if(c){if(c.tagName)return x(c);if(y(c,"array")&&a.set)return a.set.apply(a,c);if(c instanceof e)return c;if(null==t)return c=G.doc.querySelector(c),x(c)}return new s(null==c?"100%":c,null==t?"100%":t)}function v(c,a){if(a){"#text"==c&&(c=G.doc.createTextNode(a.text||""));"string"==typeof c&&(c=v(c));if("string"==typeof a)return"xlink:"==a.substring(0,6)?c.getAttributeNS(m,a.substring(6)):"xml:"==a.substring(0,4)?c.getAttributeNS(la,
a.substring(4)):c.getAttribute(a);for(var da in a)if(a[h](da)){var b=J(a[da]);b?"xlink:"==da.substring(0,6)?c.setAttributeNS(m,da.substring(6),b):"xml:"==da.substring(0,4)?c.setAttributeNS(la,da.substring(4),b):c.setAttribute(da,b):c.removeAttribute(da)}}else c=G.doc.createElementNS(la,c);return c}function y(c,a){a=J.prototype.toLowerCase.call(a);return"finite"==a?isFinite(c):"array"==a&&(c instanceof Array||Array.isArray&&Array.isArray(c))?!0:"null"==a&&null===c||a==typeof c&&null!==c||"object"==
a&&c===Object(c)||$.call(c).slice(8,-1).toLowerCase()==a}function M(c){if("function"==typeof c||Object(c)!==c)return c;var a=new c.constructor,b;for(b in c)c[h](b)&&(a[b]=M(c[b]));return a}function A(c,a,b){function m(){var e=Array.prototype.slice.call(arguments,0),f=e.join("\u2400"),d=m.cache=m.cache||{},l=m.count=m.count||[];if(d[h](f)){a:for(var e=l,l=f,B=0,H=e.length;B<H;B++)if(e[B]===l){e.push(e.splice(B,1)[0]);break a}return b?b(d[f]):d[f]}1E3<=l.length&&delete d[l.shift()];l.push(f);d[f]=c.apply(a,
e);return b?b(d[f]):d[f]}return m}function w(c,a,b,m,e,f){return null==e?(c-=b,a-=m,c||a?(180*I.atan2(-a,-c)/C+540)%360:0):w(c,a,e,f)-w(b,m,e,f)}function z(c){return c%360*C/180}function d(c){var a=[];c=c.replace(/(?:^|\s)(\w+)\(([^)]+)\)/g,function(c,b,m){m=m.split(/\s*,\s*|\s+/);"rotate"==b&&1==m.length&&m.push(0,0);"scale"==b&&(2<m.length?m=m.slice(0,2):2==m.length&&m.push(0,0),1==m.length&&m.push(m[0],0,0));"skewX"==b?a.push(["m",1,0,I.tan(z(m[0])),1,0,0]):"skewY"==b?a.push(["m",1,I.tan(z(m[0])),
0,1,0,0]):a.push([b.charAt(0)].concat(m));return c});return a}function f(c,t){var b=O(c),m=new a.Matrix;if(b)for(var e=0,f=b.length;e<f;e++){var h=b[e],d=h.length,B=J(h[0]).toLowerCase(),H=h[0]!=B,l=H?m.invert():0,E;"t"==B&&2==d?m.translate(h[1],0):"t"==B&&3==d?H?(d=l.x(0,0),B=l.y(0,0),H=l.x(h[1],h[2]),l=l.y(h[1],h[2]),m.translate(H-d,l-B)):m.translate(h[1],h[2]):"r"==B?2==d?(E=E||t,m.rotate(h[1],E.x+E.width/2,E.y+E.height/2)):4==d&&(H?(H=l.x(h[2],h[3]),l=l.y(h[2],h[3]),m.rotate(h[1],H,l)):m.rotate(h[1],
h[2],h[3])):"s"==B?2==d||3==d?(E=E||t,m.scale(h[1],h[d-1],E.x+E.width/2,E.y+E.height/2)):4==d?H?(H=l.x(h[2],h[3]),l=l.y(h[2],h[3]),m.scale(h[1],h[1],H,l)):m.scale(h[1],h[1],h[2],h[3]):5==d&&(H?(H=l.x(h[3],h[4]),l=l.y(h[3],h[4]),m.scale(h[1],h[2],H,l)):m.scale(h[1],h[2],h[3],h[4])):"m"==B&&7==d&&m.add(h[1],h[2],h[3],h[4],h[5],h[6])}return m}function n(c,t){if(null==t){var m=!0;t="linearGradient"==c.type||"radialGradient"==c.type?c.node.getAttribute("gradientTransform"):"pattern"==c.type?c.node.getAttribute("patternTransform"):
c.node.getAttribute("transform");if(!t)return new a.Matrix;t=d(t)}else t=a._.rgTransform.test(t)?J(t).replace(/\.{3}|\u2026/g,c._.transform||aa):d(t),y(t,"array")&&(t=a.path?a.path.toString.call(t):J(t)),c._.transform=t;var b=f(t,c.getBBox(1));if(m)return b;c.matrix=b}function u(c){c=c.node.ownerSVGElement&&x(c.node.ownerSVGElement)||c.node.parentNode&&x(c.node.parentNode)||a.select("svg")||a(0,0);var t=c.select("defs"),t=null==t?!1:t.node;t||(t=r("defs",c.node).node);return t}function p(c){return c.node.ownerSVGElement&&
x(c.node.ownerSVGElement)||a.select("svg")}function b(c,a,m){function b(c){if(null==c)return aa;if(c==+c)return c;v(B,{width:c});try{return B.getBBox().width}catch(a){return 0}}function h(c){if(null==c)return aa;if(c==+c)return c;v(B,{height:c});try{return B.getBBox().height}catch(a){return 0}}function e(b,B){null==a?d[b]=B(c.attr(b)||0):b==a&&(d=B(null==m?c.attr(b)||0:m))}var f=p(c).node,d={},B=f.querySelector(".svg---mgr");B||(B=v("rect"),v(B,{x:-9E9,y:-9E9,width:10,height:10,"class":"svg---mgr",
fill:"none"}),f.appendChild(B));switch(c.type){case "rect":e("rx",b),e("ry",h);case "image":e("width",b),e("height",h);case "text":e("x",b);e("y",h);break;case "circle":e("cx",b);e("cy",h);e("r",b);break;case "ellipse":e("cx",b);e("cy",h);e("rx",b);e("ry",h);break;case "line":e("x1",b);e("x2",b);e("y1",h);e("y2",h);break;case "marker":e("refX",b);e("markerWidth",b);e("refY",h);e("markerHeight",h);break;case "radialGradient":e("fx",b);e("fy",h);break;case "tspan":e("dx",b);e("dy",h);break;default:e(a,
b)}f.removeChild(B);return d}function q(c){y(c,"array")||(c=Array.prototype.slice.call(arguments,0));for(var a=0,b=0,m=this.node;this[a];)delete this[a++];for(a=0;a<c.length;a++)"set"==c[a].type?c[a].forEach(function(c){m.appendChild(c.node)}):m.appendChild(c[a].node);for(var h=m.childNodes,a=0;a<h.length;a++)this[b++]=x(h[a]);return this}function e(c){if(c.snap in E)return E[c.snap];var a=this.id=V(),b;try{b=c.ownerSVGElement}catch(m){}this.node=c;b&&(this.paper=new s(b));this.type=c.tagName;this.anims=
{};this._={transform:[]};c.snap=a;E[a]=this;"g"==this.type&&(this.add=q);if(this.type in{g:1,mask:1,pattern:1})for(var e in s.prototype)s.prototype[h](e)&&(this[e]=s.prototype[e])}function l(c){this.node=c}function r(c,a){var b=v(c);a.appendChild(b);return x(b)}function s(c,a){var b,m,f,d=s.prototype;if(c&&"svg"==c.tagName){if(c.snap in E)return E[c.snap];var l=c.ownerDocument;b=new e(c);m=c.getElementsByTagName("desc")[0];f=c.getElementsByTagName("defs")[0];m||(m=v("desc"),m.appendChild(l.createTextNode("Created with Snap")),
b.node.appendChild(m));f||(f=v("defs"),b.node.appendChild(f));b.defs=f;for(var ca in d)d[h](ca)&&(b[ca]=d[ca]);b.paper=b.root=b}else b=r("svg",G.doc.body),v(b.node,{height:a,version:1.1,width:c,xmlns:la});return b}function x(c){return!c||c instanceof e||c instanceof l?c:c.tagName&&"svg"==c.tagName.toLowerCase()?new s(c):c.tagName&&"object"==c.tagName.toLowerCase()&&"image/svg+xml"==c.type?new s(c.contentDocument.getElementsByTagName("svg")[0]):new e(c)}a.version="0.3.0";a.toString=function(){return"Snap v"+
this.version};a._={};var G={win:N,doc:N.document};a._.glob=G;var h="hasOwnProperty",J=String,K=parseFloat,U=parseInt,I=Math,P=I.max,Q=I.min,Y=I.abs,C=I.PI,aa="",$=Object.prototype.toString,F=/^\s*((#[a-f\d]{6})|(#[a-f\d]{3})|rgba?\(\s*([\d\.]+%?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+%?(?:\s*,\s*[\d\.]+%?)?)\s*\)|hsba?\(\s*([\d\.]+(?:deg|\xb0|%)?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+(?:%?\s*,\s*[\d\.]+)?%?)\s*\)|hsla?\(\s*([\d\.]+(?:deg|\xb0|%)?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+(?:%?\s*,\s*[\d\.]+)?%?)\s*\))\s*$/i;a._.separator=
RegExp("[,\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]+");var S=RegExp("[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*"),X={hs:1,rg:1},W=RegExp("([a-z])[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029,]*((-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*)+)",
"ig"),ma=RegExp("([rstm])[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029,]*((-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*)+)","ig"),Z=RegExp("(-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?)[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*",
"ig"),na=0,ba="S"+(+new Date).toString(36),V=function(){return ba+(na++).toString(36)},m="http://www.w3.org/1999/xlink",la="http://www.w3.org/2000/svg",E={},ca=a.url=function(c){return"url('#"+c+"')"};a._.$=v;a._.id=V;a.format=function(){var c=/\{([^\}]+)\}/g,a=/(?:(?:^|\.)(.+?)(?=\[|\.|$|\()|\[('|")(.+?)\2\])(\(\))?/g,b=function(c,b,m){var h=m;b.replace(a,function(c,a,b,m,t){a=a||m;h&&(a in h&&(h=h[a]),"function"==typeof h&&t&&(h=h()))});return h=(null==h||h==m?c:h)+""};return function(a,m){return J(a).replace(c,
function(c,a){return b(c,a,m)})}}();a._.clone=M;a._.cacher=A;a.rad=z;a.deg=function(c){return 180*c/C%360};a.angle=w;a.is=y;a.snapTo=function(c,a,b){b=y(b,"finite")?b:10;if(y(c,"array"))for(var m=c.length;m--;){if(Y(c[m]-a)<=b)return c[m]}else{c=+c;m=a%c;if(m<b)return a-m;if(m>c-b)return a-m+c}return a};a.getRGB=A(function(c){if(!c||(c=J(c)).indexOf("-")+1)return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka};if("none"==c)return{r:-1,g:-1,b:-1,hex:"none",toString:ka};!X[h](c.toLowerCase().substring(0,
2))&&"#"!=c.charAt()&&(c=T(c));if(!c)return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka};var b,m,e,f,d;if(c=c.match(F)){c[2]&&(e=U(c[2].substring(5),16),m=U(c[2].substring(3,5),16),b=U(c[2].substring(1,3),16));c[3]&&(e=U((d=c[3].charAt(3))+d,16),m=U((d=c[3].charAt(2))+d,16),b=U((d=c[3].charAt(1))+d,16));c[4]&&(d=c[4].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b*=2.55),m=K(d[1]),"%"==d[1].slice(-1)&&(m*=2.55),e=K(d[2]),"%"==d[2].slice(-1)&&(e*=2.55),"rgba"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),
d[3]&&"%"==d[3].slice(-1)&&(f/=100));if(c[5])return d=c[5].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b/=100),m=K(d[1]),"%"==d[1].slice(-1)&&(m/=100),e=K(d[2]),"%"==d[2].slice(-1)&&(e/=100),"deg"!=d[0].slice(-3)&&"\u00b0"!=d[0].slice(-1)||(b/=360),"hsba"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),d[3]&&"%"==d[3].slice(-1)&&(f/=100),a.hsb2rgb(b,m,e,f);if(c[6])return d=c[6].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b/=100),m=K(d[1]),"%"==d[1].slice(-1)&&(m/=100),e=K(d[2]),"%"==d[2].slice(-1)&&(e/=100),
"deg"!=d[0].slice(-3)&&"\u00b0"!=d[0].slice(-1)||(b/=360),"hsla"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),d[3]&&"%"==d[3].slice(-1)&&(f/=100),a.hsl2rgb(b,m,e,f);b=Q(I.round(b),255);m=Q(I.round(m),255);e=Q(I.round(e),255);f=Q(P(f,0),1);c={r:b,g:m,b:e,toString:ka};c.hex="#"+(16777216|e|m<<8|b<<16).toString(16).slice(1);c.opacity=y(f,"finite")?f:1;return c}return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka}},a);a.hsb=A(function(c,b,m){return a.hsb2rgb(c,b,m).hex});a.hsl=A(function(c,b,m){return a.hsl2rgb(c,
b,m).hex});a.rgb=A(function(c,a,b,m){if(y(m,"finite")){var e=I.round;return"rgba("+[e(c),e(a),e(b),+m.toFixed(2)]+")"}return"#"+(16777216|b|a<<8|c<<16).toString(16).slice(1)});var T=function(c){var a=G.doc.getElementsByTagName("head")[0]||G.doc.getElementsByTagName("svg")[0];T=A(function(c){if("red"==c.toLowerCase())return"rgb(255, 0, 0)";a.style.color="rgb(255, 0, 0)";a.style.color=c;c=G.doc.defaultView.getComputedStyle(a,aa).getPropertyValue("color");return"rgb(255, 0, 0)"==c?null:c});return T(c)},
qa=function(){return"hsb("+[this.h,this.s,this.b]+")"},ra=function(){return"hsl("+[this.h,this.s,this.l]+")"},ka=function(){return 1==this.opacity||null==this.opacity?this.hex:"rgba("+[this.r,this.g,this.b,this.opacity]+")"},D=function(c,b,m){null==b&&y(c,"object")&&"r"in c&&"g"in c&&"b"in c&&(m=c.b,b=c.g,c=c.r);null==b&&y(c,string)&&(m=a.getRGB(c),c=m.r,b=m.g,m=m.b);if(1<c||1<b||1<m)c/=255,b/=255,m/=255;return[c,b,m]},oa=function(c,b,m,e){c=I.round(255*c);b=I.round(255*b);m=I.round(255*m);c={r:c,
g:b,b:m,opacity:y(e,"finite")?e:1,hex:a.rgb(c,b,m),toString:ka};y(e,"finite")&&(c.opacity=e);return c};a.color=function(c){var b;y(c,"object")&&"h"in c&&"s"in c&&"b"in c?(b=a.hsb2rgb(c),c.r=b.r,c.g=b.g,c.b=b.b,c.opacity=1,c.hex=b.hex):y(c,"object")&&"h"in c&&"s"in c&&"l"in c?(b=a.hsl2rgb(c),c.r=b.r,c.g=b.g,c.b=b.b,c.opacity=1,c.hex=b.hex):(y(c,"string")&&(c=a.getRGB(c)),y(c,"object")&&"r"in c&&"g"in c&&"b"in c&&!("error"in c)?(b=a.rgb2hsl(c),c.h=b.h,c.s=b.s,c.l=b.l,b=a.rgb2hsb(c),c.v=b.b):(c={hex:"none"},
c.r=c.g=c.b=c.h=c.s=c.v=c.l=-1,c.error=1));c.toString=ka;return c};a.hsb2rgb=function(c,a,b,m){y(c,"object")&&"h"in c&&"s"in c&&"b"in c&&(b=c.b,a=c.s,c=c.h,m=c.o);var e,h,d;c=360*c%360/60;d=b*a;a=d*(1-Y(c%2-1));b=e=h=b-d;c=~~c;b+=[d,a,0,0,a,d][c];e+=[a,d,d,a,0,0][c];h+=[0,0,a,d,d,a][c];return oa(b,e,h,m)};a.hsl2rgb=function(c,a,b,m){y(c,"object")&&"h"in c&&"s"in c&&"l"in c&&(b=c.l,a=c.s,c=c.h);if(1<c||1<a||1<b)c/=360,a/=100,b/=100;var e,h,d;c=360*c%360/60;d=2*a*(0.5>b?b:1-b);a=d*(1-Y(c%2-1));b=e=
h=b-d/2;c=~~c;b+=[d,a,0,0,a,d][c];e+=[a,d,d,a,0,0][c];h+=[0,0,a,d,d,a][c];return oa(b,e,h,m)};a.rgb2hsb=function(c,a,b){b=D(c,a,b);c=b[0];a=b[1];b=b[2];var m,e;m=P(c,a,b);e=m-Q(c,a,b);c=((0==e?0:m==c?(a-b)/e:m==a?(b-c)/e+2:(c-a)/e+4)+360)%6*60/360;return{h:c,s:0==e?0:e/m,b:m,toString:qa}};a.rgb2hsl=function(c,a,b){b=D(c,a,b);c=b[0];a=b[1];b=b[2];var m,e,h;m=P(c,a,b);e=Q(c,a,b);h=m-e;c=((0==h?0:m==c?(a-b)/h:m==a?(b-c)/h+2:(c-a)/h+4)+360)%6*60/360;m=(m+e)/2;return{h:c,s:0==h?0:0.5>m?h/(2*m):h/(2-2*
m),l:m,toString:ra}};a.parsePathString=function(c){if(!c)return null;var b=a.path(c);if(b.arr)return a.path.clone(b.arr);var m={a:7,c:6,o:2,h:1,l:2,m:2,r:4,q:4,s:4,t:2,v:1,u:3,z:0},e=[];y(c,"array")&&y(c[0],"array")&&(e=a.path.clone(c));e.length||J(c).replace(W,function(c,a,b){var h=[];c=a.toLowerCase();b.replace(Z,function(c,a){a&&h.push(+a)});"m"==c&&2<h.length&&(e.push([a].concat(h.splice(0,2))),c="l",a="m"==a?"l":"L");"o"==c&&1==h.length&&e.push([a,h[0] ]);if("r"==c)e.push([a].concat(h));else for(;h.length>=
m[c]&&(e.push([a].concat(h.splice(0,m[c]))),m[c]););});e.toString=a.path.toString;b.arr=a.path.clone(e);return e};var O=a.parseTransformString=function(c){if(!c)return null;var b=[];y(c,"array")&&y(c[0],"array")&&(b=a.path.clone(c));b.length||J(c).replace(ma,function(c,a,m){var e=[];a.toLowerCase();m.replace(Z,function(c,a){a&&e.push(+a)});b.push([a].concat(e))});b.toString=a.path.toString;return b};a._.svgTransform2string=d;a._.rgTransform=RegExp("^[a-z][\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*-?\\.?\\d",
"i");a._.transform2matrix=f;a._unit2px=b;a._.getSomeDefs=u;a._.getSomeSVG=p;a.select=function(c){return x(G.doc.querySelector(c))};a.selectAll=function(c){c=G.doc.querySelectorAll(c);for(var b=(a.set||Array)(),m=0;m<c.length;m++)b.push(x(c[m]));return b};setInterval(function(){for(var c in E)if(E[h](c)){var a=E[c],b=a.node;("svg"!=a.type&&!b.ownerSVGElement||"svg"==a.type&&(!b.parentNode||"ownerSVGElement"in b.parentNode&&!b.ownerSVGElement))&&delete E[c]}},1E4);(function(c){function m(c){function a(c,
b){var m=v(c.node,b);(m=(m=m&&m.match(d))&&m[2])&&"#"==m.charAt()&&(m=m.substring(1))&&(f[m]=(f[m]||[]).concat(function(a){var m={};m[b]=ca(a);v(c.node,m)}))}function b(c){var a=v(c.node,"xlink:href");a&&"#"==a.charAt()&&(a=a.substring(1))&&(f[a]=(f[a]||[]).concat(function(a){c.attr("xlink:href","#"+a)}))}var e=c.selectAll("*"),h,d=/^\s*url\(("|'|)(.*)\1\)\s*$/;c=[];for(var f={},l=0,E=e.length;l<E;l++){h=e[l];a(h,"fill");a(h,"stroke");a(h,"filter");a(h,"mask");a(h,"clip-path");b(h);var t=v(h.node,
"id");t&&(v(h.node,{id:h.id}),c.push({old:t,id:h.id}))}l=0;for(E=c.length;l<E;l++)if(e=f[c[l].old])for(h=0,t=e.length;h<t;h++)e[h](c[l].id)}function e(c,a,b){return function(m){m=m.slice(c,a);1==m.length&&(m=m[0]);return b?b(m):m}}function d(c){return function(){var a=c?"<"+this.type:"",b=this.node.attributes,m=this.node.childNodes;if(c)for(var e=0,h=b.length;e<h;e++)a+=" "+b[e].name+'="'+b[e].value.replace(/"/g,'\\"')+'"';if(m.length){c&&(a+=">");e=0;for(h=m.length;e<h;e++)3==m[e].nodeType?a+=m[e].nodeValue:
1==m[e].nodeType&&(a+=x(m[e]).toString());c&&(a+="</"+this.type+">")}else c&&(a+="/>");return a}}c.attr=function(c,a){if(!c)return this;if(y(c,"string"))if(1<arguments.length){var b={};b[c]=a;c=b}else return k("snap.util.getattr."+c,this).firstDefined();for(var m in c)c[h](m)&&k("snap.util.attr."+m,this,c[m]);return this};c.getBBox=function(c){if(!a.Matrix||!a.path)return this.node.getBBox();var b=this,m=new a.Matrix;if(b.removed)return a._.box();for(;"use"==b.type;)if(c||(m=m.add(b.transform().localMatrix.translate(b.attr("x")||
0,b.attr("y")||0))),b.original)b=b.original;else var e=b.attr("xlink:href"),b=b.original=b.node.ownerDocument.getElementById(e.substring(e.indexOf("#")+1));var e=b._,h=a.path.get[b.type]||a.path.get.deflt;try{if(c)return e.bboxwt=h?a.path.getBBox(b.realPath=h(b)):a._.box(b.node.getBBox()),a._.box(e.bboxwt);b.realPath=h(b);b.matrix=b.transform().localMatrix;e.bbox=a.path.getBBox(a.path.map(b.realPath,m.add(b.matrix)));return a._.box(e.bbox)}catch(d){return a._.box()}};var f=function(){return this.string};
c.transform=function(c){var b=this._;if(null==c){var m=this;c=new a.Matrix(this.node.getCTM());for(var e=n(this),h=[e],d=new a.Matrix,l=e.toTransformString(),b=J(e)==J(this.matrix)?J(b.transform):l;"svg"!=m.type&&(m=m.parent());)h.push(n(m));for(m=h.length;m--;)d.add(h[m]);return{string:b,globalMatrix:c,totalMatrix:d,localMatrix:e,diffMatrix:c.clone().add(e.invert()),global:c.toTransformString(),total:d.toTransformString(),local:l,toString:f}}c instanceof a.Matrix?this.matrix=c:n(this,c);this.node&&
("linearGradient"==this.type||"radialGradient"==this.type?v(this.node,{gradientTransform:this.matrix}):"pattern"==this.type?v(this.node,{patternTransform:this.matrix}):v(this.node,{transform:this.matrix}));return this};c.parent=function(){return x(this.node.parentNode)};c.append=c.add=function(c){if(c){if("set"==c.type){var a=this;c.forEach(function(c){a.add(c)});return this}c=x(c);this.node.appendChild(c.node);c.paper=this.paper}return this};c.appendTo=function(c){c&&(c=x(c),c.append(this));return this};
c.prepend=function(c){if(c){if("set"==c.type){var a=this,b;c.forEach(function(c){b?b.after(c):a.prepend(c);b=c});return this}c=x(c);var m=c.parent();this.node.insertBefore(c.node,this.node.firstChild);this.add&&this.add();c.paper=this.paper;this.parent()&&this.parent().add();m&&m.add()}return this};c.prependTo=function(c){c=x(c);c.prepend(this);return this};c.before=function(c){if("set"==c.type){var a=this;c.forEach(function(c){var b=c.parent();a.node.parentNode.insertBefore(c.node,a.node);b&&b.add()});
this.parent().add();return this}c=x(c);var b=c.parent();this.node.parentNode.insertBefore(c.node,this.node);this.parent()&&this.parent().add();b&&b.add();c.paper=this.paper;return this};c.after=function(c){c=x(c);var a=c.parent();this.node.nextSibling?this.node.parentNode.insertBefore(c.node,this.node.nextSibling):this.node.parentNode.appendChild(c.node);this.parent()&&this.parent().add();a&&a.add();c.paper=this.paper;return this};c.insertBefore=function(c){c=x(c);var a=this.parent();c.node.parentNode.insertBefore(this.node,
c.node);this.paper=c.paper;a&&a.add();c.parent()&&c.parent().add();return this};c.insertAfter=function(c){c=x(c);var a=this.parent();c.node.parentNode.insertBefore(this.node,c.node.nextSibling);this.paper=c.paper;a&&a.add();c.parent()&&c.parent().add();return this};c.remove=function(){var c=this.parent();this.node.parentNode&&this.node.parentNode.removeChild(this.node);delete this.paper;this.removed=!0;c&&c.add();return this};c.select=function(c){return x(this.node.querySelector(c))};c.selectAll=
function(c){c=this.node.querySelectorAll(c);for(var b=(a.set||Array)(),m=0;m<c.length;m++)b.push(x(c[m]));return b};c.asPX=function(c,a){null==a&&(a=this.attr(c));return+b(this,c,a)};c.use=function(){var c,a=this.node.id;a||(a=this.id,v(this.node,{id:a}));c="linearGradient"==this.type||"radialGradient"==this.type||"pattern"==this.type?r(this.type,this.node.parentNode):r("use",this.node.parentNode);v(c.node,{"xlink:href":"#"+a});c.original=this;return c};var l=/\S+/g;c.addClass=function(c){var a=(c||
"").match(l)||[];c=this.node;var b=c.className.baseVal,m=b.match(l)||[],e,h,d;if(a.length){for(e=0;d=a[e++];)h=m.indexOf(d),~h||m.push(d);a=m.join(" ");b!=a&&(c.className.baseVal=a)}return this};c.removeClass=function(c){var a=(c||"").match(l)||[];c=this.node;var b=c.className.baseVal,m=b.match(l)||[],e,h;if(m.length){for(e=0;h=a[e++];)h=m.indexOf(h),~h&&m.splice(h,1);a=m.join(" ");b!=a&&(c.className.baseVal=a)}return this};c.hasClass=function(c){return!!~(this.node.className.baseVal.match(l)||[]).indexOf(c)};
c.toggleClass=function(c,a){if(null!=a)return a?this.addClass(c):this.removeClass(c);var b=(c||"").match(l)||[],m=this.node,e=m.className.baseVal,h=e.match(l)||[],d,f,E;for(d=0;E=b[d++];)f=h.indexOf(E),~f?h.splice(f,1):h.push(E);b=h.join(" ");e!=b&&(m.className.baseVal=b);return this};c.clone=function(){var c=x(this.node.cloneNode(!0));v(c.node,"id")&&v(c.node,{id:c.id});m(c);c.insertAfter(this);return c};c.toDefs=function(){u(this).appendChild(this.node);return this};c.pattern=c.toPattern=function(c,
a,b,m){var e=r("pattern",u(this));null==c&&(c=this.getBBox());y(c,"object")&&"x"in c&&(a=c.y,b=c.width,m=c.height,c=c.x);v(e.node,{x:c,y:a,width:b,height:m,patternUnits:"userSpaceOnUse",id:e.id,viewBox:[c,a,b,m].join(" ")});e.node.appendChild(this.node);return e};c.marker=function(c,a,b,m,e,h){var d=r("marker",u(this));null==c&&(c=this.getBBox());y(c,"object")&&"x"in c&&(a=c.y,b=c.width,m=c.height,e=c.refX||c.cx,h=c.refY||c.cy,c=c.x);v(d.node,{viewBox:[c,a,b,m].join(" "),markerWidth:b,markerHeight:m,
orient:"auto",refX:e||0,refY:h||0,id:d.id});d.node.appendChild(this.node);return d};var E=function(c,a,b,m){"function"!=typeof b||b.length||(m=b,b=L.linear);this.attr=c;this.dur=a;b&&(this.easing=b);m&&(this.callback=m)};a._.Animation=E;a.animation=function(c,a,b,m){return new E(c,a,b,m)};c.inAnim=function(){var c=[],a;for(a in this.anims)this.anims[h](a)&&function(a){c.push({anim:new E(a._attrs,a.dur,a.easing,a._callback),mina:a,curStatus:a.status(),status:function(c){return a.status(c)},stop:function(){a.stop()}})}(this.anims[a]);
return c};a.animate=function(c,a,b,m,e,h){"function"!=typeof e||e.length||(h=e,e=L.linear);var d=L.time();c=L(c,a,d,d+m,L.time,b,e);h&&k.once("mina.finish."+c.id,h);return c};c.stop=function(){for(var c=this.inAnim(),a=0,b=c.length;a<b;a++)c[a].stop();return this};c.animate=function(c,a,b,m){"function"!=typeof b||b.length||(m=b,b=L.linear);c instanceof E&&(m=c.callback,b=c.easing,a=b.dur,c=c.attr);var d=[],f=[],l={},t,ca,n,T=this,q;for(q in c)if(c[h](q)){T.equal?(n=T.equal(q,J(c[q])),t=n.from,ca=
n.to,n=n.f):(t=+T.attr(q),ca=+c[q]);var la=y(t,"array")?t.length:1;l[q]=e(d.length,d.length+la,n);d=d.concat(t);f=f.concat(ca)}t=L.time();var p=L(d,f,t,t+a,L.time,function(c){var a={},b;for(b in l)l[h](b)&&(a[b]=l[b](c));T.attr(a)},b);T.anims[p.id]=p;p._attrs=c;p._callback=m;k("snap.animcreated."+T.id,p);k.once("mina.finish."+p.id,function(){delete T.anims[p.id];m&&m.call(T)});k.once("mina.stop."+p.id,function(){delete T.anims[p.id]});return T};var T={};c.data=function(c,b){var m=T[this.id]=T[this.id]||
{};if(0==arguments.length)return k("snap.data.get."+this.id,this,m,null),m;if(1==arguments.length){if(a.is(c,"object")){for(var e in c)c[h](e)&&this.data(e,c[e]);return this}k("snap.data.get."+this.id,this,m[c],c);return m[c]}m[c]=b;k("snap.data.set."+this.id,this,b,c);return this};c.removeData=function(c){null==c?T[this.id]={}:T[this.id]&&delete T[this.id][c];return this};c.outerSVG=c.toString=d(1);c.innerSVG=d()})(e.prototype);a.parse=function(c){var a=G.doc.createDocumentFragment(),b=!0,m=G.doc.createElement("div");
c=J(c);c.match(/^\s*<\s*svg(?:\s|>)/)||(c="<svg>"+c+"</svg>",b=!1);m.innerHTML=c;if(c=m.getElementsByTagName("svg")[0])if(b)a=c;else for(;c.firstChild;)a.appendChild(c.firstChild);m.innerHTML=aa;return new l(a)};l.prototype.select=e.prototype.select;l.prototype.selectAll=e.prototype.selectAll;a.fragment=function(){for(var c=Array.prototype.slice.call(arguments,0),b=G.doc.createDocumentFragment(),m=0,e=c.length;m<e;m++){var h=c[m];h.node&&h.node.nodeType&&b.appendChild(h.node);h.nodeType&&b.appendChild(h);
"string"==typeof h&&b.appendChild(a.parse(h).node)}return new l(b)};a._.make=r;a._.wrap=x;s.prototype.el=function(c,a){var b=r(c,this.node);a&&b.attr(a);return b};k.on("snap.util.getattr",function(){var c=k.nt(),c=c.substring(c.lastIndexOf(".")+1),a=c.replace(/[A-Z]/g,function(c){return"-"+c.toLowerCase()});return pa[h](a)?this.node.ownerDocument.defaultView.getComputedStyle(this.node,null).getPropertyValue(a):v(this.node,c)});var pa={"alignment-baseline":0,"baseline-shift":0,clip:0,"clip-path":0,
"clip-rule":0,color:0,"color-interpolation":0,"color-interpolation-filters":0,"color-profile":0,"color-rendering":0,cursor:0,direction:0,display:0,"dominant-baseline":0,"enable-background":0,fill:0,"fill-opacity":0,"fill-rule":0,filter:0,"flood-color":0,"flood-opacity":0,font:0,"font-family":0,"font-size":0,"font-size-adjust":0,"font-stretch":0,"font-style":0,"font-variant":0,"font-weight":0,"glyph-orientation-horizontal":0,"glyph-orientation-vertical":0,"image-rendering":0,kerning:0,"letter-spacing":0,
"lighting-color":0,marker:0,"marker-end":0,"marker-mid":0,"marker-start":0,mask:0,opacity:0,overflow:0,"pointer-events":0,"shape-rendering":0,"stop-color":0,"stop-opacity":0,stroke:0,"stroke-dasharray":0,"stroke-dashoffset":0,"stroke-linecap":0,"stroke-linejoin":0,"stroke-miterlimit":0,"stroke-opacity":0,"stroke-width":0,"text-anchor":0,"text-decoration":0,"text-rendering":0,"unicode-bidi":0,visibility:0,"word-spacing":0,"writing-mode":0};k.on("snap.util.attr",function(c){var a=k.nt(),b={},a=a.substring(a.lastIndexOf(".")+
1);b[a]=c;var m=a.replace(/-(\w)/gi,function(c,a){return a.toUpperCase()}),a=a.replace(/[A-Z]/g,function(c){return"-"+c.toLowerCase()});pa[h](a)?this.node.style[m]=null==c?aa:c:v(this.node,b)});a.ajax=function(c,a,b,m){var e=new XMLHttpRequest,h=V();if(e){if(y(a,"function"))m=b,b=a,a=null;else if(y(a,"object")){var d=[],f;for(f in a)a.hasOwnProperty(f)&&d.push(encodeURIComponent(f)+"="+encodeURIComponent(a[f]));a=d.join("&")}e.open(a?"POST":"GET",c,!0);a&&(e.setRequestHeader("X-Requested-With","XMLHttpRequest"),
e.setRequestHeader("Content-type","application/x-www-form-urlencoded"));b&&(k.once("snap.ajax."+h+".0",b),k.once("snap.ajax."+h+".200",b),k.once("snap.ajax."+h+".304",b));e.onreadystatechange=function(){4==e.readyState&&k("snap.ajax."+h+"."+e.status,m,e)};if(4==e.readyState)return e;e.send(a);return e}};a.load=function(c,b,m){a.ajax(c,function(c){c=a.parse(c.responseText);m?b.call(m,c):b(c)})};a.getElementByPoint=function(c,a){var b,m,e=G.doc.elementFromPoint(c,a);if(G.win.opera&&"svg"==e.tagName){b=
e;m=b.getBoundingClientRect();b=b.ownerDocument;var h=b.body,d=b.documentElement;b=m.top+(g.win.pageYOffset||d.scrollTop||h.scrollTop)-(d.clientTop||h.clientTop||0);m=m.left+(g.win.pageXOffset||d.scrollLeft||h.scrollLeft)-(d.clientLeft||h.clientLeft||0);h=e.createSVGRect();h.x=c-m;h.y=a-b;h.width=h.height=1;b=e.getIntersectionList(h,null);b.length&&(e=b[b.length-1])}return e?x(e):null};a.plugin=function(c){c(a,e,s,G,l)};return G.win.Snap=a}();C.plugin(function(a,k,y,M,A){function w(a,d,f,b,q,e){null==
d&&"[object SVGMatrix]"==z.call(a)?(this.a=a.a,this.b=a.b,this.c=a.c,this.d=a.d,this.e=a.e,this.f=a.f):null!=a?(this.a=+a,this.b=+d,this.c=+f,this.d=+b,this.e=+q,this.f=+e):(this.a=1,this.c=this.b=0,this.d=1,this.f=this.e=0)}var z=Object.prototype.toString,d=String,f=Math;(function(n){function k(a){return a[0]*a[0]+a[1]*a[1]}function p(a){var d=f.sqrt(k(a));a[0]&&(a[0]/=d);a[1]&&(a[1]/=d)}n.add=function(a,d,e,f,n,p){var k=[[],[],[] ],u=[[this.a,this.c,this.e],[this.b,this.d,this.f],[0,0,1] ];d=[[a,
e,n],[d,f,p],[0,0,1] ];a&&a instanceof w&&(d=[[a.a,a.c,a.e],[a.b,a.d,a.f],[0,0,1] ]);for(a=0;3>a;a++)for(e=0;3>e;e++){for(f=n=0;3>f;f++)n+=u[a][f]*d[f][e];k[a][e]=n}this.a=k[0][0];this.b=k[1][0];this.c=k[0][1];this.d=k[1][1];this.e=k[0][2];this.f=k[1][2];return this};n.invert=function(){var a=this.a*this.d-this.b*this.c;return new w(this.d/a,-this.b/a,-this.c/a,this.a/a,(this.c*this.f-this.d*this.e)/a,(this.b*this.e-this.a*this.f)/a)};n.clone=function(){return new w(this.a,this.b,this.c,this.d,this.e,
this.f)};n.translate=function(a,d){return this.add(1,0,0,1,a,d)};n.scale=function(a,d,e,f){null==d&&(d=a);(e||f)&&this.add(1,0,0,1,e,f);this.add(a,0,0,d,0,0);(e||f)&&this.add(1,0,0,1,-e,-f);return this};n.rotate=function(b,d,e){b=a.rad(b);d=d||0;e=e||0;var l=+f.cos(b).toFixed(9);b=+f.sin(b).toFixed(9);this.add(l,b,-b,l,d,e);return this.add(1,0,0,1,-d,-e)};n.x=function(a,d){return a*this.a+d*this.c+this.e};n.y=function(a,d){return a*this.b+d*this.d+this.f};n.get=function(a){return+this[d.fromCharCode(97+
a)].toFixed(4)};n.toString=function(){return"matrix("+[this.get(0),this.get(1),this.get(2),this.get(3),this.get(4),this.get(5)].join()+")"};n.offset=function(){return[this.e.toFixed(4),this.f.toFixed(4)]};n.determinant=function(){return this.a*this.d-this.b*this.c};n.split=function(){var b={};b.dx=this.e;b.dy=this.f;var d=[[this.a,this.c],[this.b,this.d] ];b.scalex=f.sqrt(k(d[0]));p(d[0]);b.shear=d[0][0]*d[1][0]+d[0][1]*d[1][1];d[1]=[d[1][0]-d[0][0]*b.shear,d[1][1]-d[0][1]*b.shear];b.scaley=f.sqrt(k(d[1]));
p(d[1]);b.shear/=b.scaley;0>this.determinant()&&(b.scalex=-b.scalex);var e=-d[0][1],d=d[1][1];0>d?(b.rotate=a.deg(f.acos(d)),0>e&&(b.rotate=360-b.rotate)):b.rotate=a.deg(f.asin(e));b.isSimple=!+b.shear.toFixed(9)&&(b.scalex.toFixed(9)==b.scaley.toFixed(9)||!b.rotate);b.isSuperSimple=!+b.shear.toFixed(9)&&b.scalex.toFixed(9)==b.scaley.toFixed(9)&&!b.rotate;b.noRotation=!+b.shear.toFixed(9)&&!b.rotate;return b};n.toTransformString=function(a){a=a||this.split();if(+a.shear.toFixed(9))return"m"+[this.get(0),
this.get(1),this.get(2),this.get(3),this.get(4),this.get(5)];a.scalex=+a.scalex.toFixed(4);a.scaley=+a.scaley.toFixed(4);a.rotate=+a.rotate.toFixed(4);return(a.dx||a.dy?"t"+[+a.dx.toFixed(4),+a.dy.toFixed(4)]:"")+(1!=a.scalex||1!=a.scaley?"s"+[a.scalex,a.scaley,0,0]:"")+(a.rotate?"r"+[+a.rotate.toFixed(4),0,0]:"")}})(w.prototype);a.Matrix=w;a.matrix=function(a,d,f,b,k,e){return new w(a,d,f,b,k,e)}});C.plugin(function(a,v,y,M,A){function w(h){return function(d){k.stop();d instanceof A&&1==d.node.childNodes.length&&
("radialGradient"==d.node.firstChild.tagName||"linearGradient"==d.node.firstChild.tagName||"pattern"==d.node.firstChild.tagName)&&(d=d.node.firstChild,b(this).appendChild(d),d=u(d));if(d instanceof v)if("radialGradient"==d.type||"linearGradient"==d.type||"pattern"==d.type){d.node.id||e(d.node,{id:d.id});var f=l(d.node.id)}else f=d.attr(h);else f=a.color(d),f.error?(f=a(b(this).ownerSVGElement).gradient(d))?(f.node.id||e(f.node,{id:f.id}),f=l(f.node.id)):f=d:f=r(f);d={};d[h]=f;e(this.node,d);this.node.style[h]=
x}}function z(a){k.stop();a==+a&&(a+="px");this.node.style.fontSize=a}function d(a){var b=[];a=a.childNodes;for(var e=0,f=a.length;e<f;e++){var l=a[e];3==l.nodeType&&b.push(l.nodeValue);"tspan"==l.tagName&&(1==l.childNodes.length&&3==l.firstChild.nodeType?b.push(l.firstChild.nodeValue):b.push(d(l)))}return b}function f(){k.stop();return this.node.style.fontSize}var n=a._.make,u=a._.wrap,p=a.is,b=a._.getSomeDefs,q=/^url\(#?([^)]+)\)$/,e=a._.$,l=a.url,r=String,s=a._.separator,x="";k.on("snap.util.attr.mask",
function(a){if(a instanceof v||a instanceof A){k.stop();a instanceof A&&1==a.node.childNodes.length&&(a=a.node.firstChild,b(this).appendChild(a),a=u(a));if("mask"==a.type)var d=a;else d=n("mask",b(this)),d.node.appendChild(a.node);!d.node.id&&e(d.node,{id:d.id});e(this.node,{mask:l(d.id)})}});(function(a){k.on("snap.util.attr.clip",a);k.on("snap.util.attr.clip-path",a);k.on("snap.util.attr.clipPath",a)})(function(a){if(a instanceof v||a instanceof A){k.stop();if("clipPath"==a.type)var d=a;else d=
n("clipPath",b(this)),d.node.appendChild(a.node),!d.node.id&&e(d.node,{id:d.id});e(this.node,{"clip-path":l(d.id)})}});k.on("snap.util.attr.fill",w("fill"));k.on("snap.util.attr.stroke",w("stroke"));var G=/^([lr])(?:\(([^)]*)\))?(.*)$/i;k.on("snap.util.grad.parse",function(a){a=r(a);var b=a.match(G);if(!b)return null;a=b[1];var e=b[2],b=b[3],e=e.split(/\s*,\s*/).map(function(a){return+a==a?+a:a});1==e.length&&0==e[0]&&(e=[]);b=b.split("-");b=b.map(function(a){a=a.split(":");var b={color:a[0]};a[1]&&
(b.offset=parseFloat(a[1]));return b});return{type:a,params:e,stops:b}});k.on("snap.util.attr.d",function(b){k.stop();p(b,"array")&&p(b[0],"array")&&(b=a.path.toString.call(b));b=r(b);b.match(/[ruo]/i)&&(b=a.path.toAbsolute(b));e(this.node,{d:b})})(-1);k.on("snap.util.attr.#text",function(a){k.stop();a=r(a);for(a=M.doc.createTextNode(a);this.node.firstChild;)this.node.removeChild(this.node.firstChild);this.node.appendChild(a)})(-1);k.on("snap.util.attr.path",function(a){k.stop();this.attr({d:a})})(-1);
k.on("snap.util.attr.class",function(a){k.stop();this.node.className.baseVal=a})(-1);k.on("snap.util.attr.viewBox",function(a){a=p(a,"object")&&"x"in a?[a.x,a.y,a.width,a.height].join(" "):p(a,"array")?a.join(" "):a;e(this.node,{viewBox:a});k.stop()})(-1);k.on("snap.util.attr.transform",function(a){this.transform(a);k.stop()})(-1);k.on("snap.util.attr.r",function(a){"rect"==this.type&&(k.stop(),e(this.node,{rx:a,ry:a}))})(-1);k.on("snap.util.attr.textpath",function(a){k.stop();if("text"==this.type){var d,
f;if(!a&&this.textPath){for(a=this.textPath;a.node.firstChild;)this.node.appendChild(a.node.firstChild);a.remove();delete this.textPath}else if(p(a,"string")?(d=b(this),a=u(d.parentNode).path(a),d.appendChild(a.node),d=a.id,a.attr({id:d})):(a=u(a),a instanceof v&&(d=a.attr("id"),d||(d=a.id,a.attr({id:d})))),d)if(a=this.textPath,f=this.node,a)a.attr({"xlink:href":"#"+d});else{for(a=e("textPath",{"xlink:href":"#"+d});f.firstChild;)a.appendChild(f.firstChild);f.appendChild(a);this.textPath=u(a)}}})(-1);
k.on("snap.util.attr.text",function(a){if("text"==this.type){for(var b=this.node,d=function(a){var b=e("tspan");if(p(a,"array"))for(var f=0;f<a.length;f++)b.appendChild(d(a[f]));else b.appendChild(M.doc.createTextNode(a));b.normalize&&b.normalize();return b};b.firstChild;)b.removeChild(b.firstChild);for(a=d(a);a.firstChild;)b.appendChild(a.firstChild)}k.stop()})(-1);k.on("snap.util.attr.fontSize",z)(-1);k.on("snap.util.attr.font-size",z)(-1);k.on("snap.util.getattr.transform",function(){k.stop();
return this.transform()})(-1);k.on("snap.util.getattr.textpath",function(){k.stop();return this.textPath})(-1);(function(){function b(d){return function(){k.stop();var b=M.doc.defaultView.getComputedStyle(this.node,null).getPropertyValue("marker-"+d);return"none"==b?b:a(M.doc.getElementById(b.match(q)[1]))}}function d(a){return function(b){k.stop();var d="marker"+a.charAt(0).toUpperCase()+a.substring(1);if(""==b||!b)this.node.style[d]="none";else if("marker"==b.type){var f=b.node.id;f||e(b.node,{id:b.id});
this.node.style[d]=l(f)}}}k.on("snap.util.getattr.marker-end",b("end"))(-1);k.on("snap.util.getattr.markerEnd",b("end"))(-1);k.on("snap.util.getattr.marker-start",b("start"))(-1);k.on("snap.util.getattr.markerStart",b("start"))(-1);k.on("snap.util.getattr.marker-mid",b("mid"))(-1);k.on("snap.util.getattr.markerMid",b("mid"))(-1);k.on("snap.util.attr.marker-end",d("end"))(-1);k.on("snap.util.attr.markerEnd",d("end"))(-1);k.on("snap.util.attr.marker-start",d("start"))(-1);k.on("snap.util.attr.markerStart",
d("start"))(-1);k.on("snap.util.attr.marker-mid",d("mid"))(-1);k.on("snap.util.attr.markerMid",d("mid"))(-1)})();k.on("snap.util.getattr.r",function(){if("rect"==this.type&&e(this.node,"rx")==e(this.node,"ry"))return k.stop(),e(this.node,"rx")})(-1);k.on("snap.util.getattr.text",function(){if("text"==this.type||"tspan"==this.type){k.stop();var a=d(this.node);return 1==a.length?a[0]:a}})(-1);k.on("snap.util.getattr.#text",function(){return this.node.textContent})(-1);k.on("snap.util.getattr.viewBox",
function(){k.stop();var b=e(this.node,"viewBox");if(b)return b=b.split(s),a._.box(+b[0],+b[1],+b[2],+b[3])})(-1);k.on("snap.util.getattr.points",function(){var a=e(this.node,"points");k.stop();if(a)return a.split(s)})(-1);k.on("snap.util.getattr.path",function(){var a=e(this.node,"d");k.stop();return a})(-1);k.on("snap.util.getattr.class",function(){return this.node.className.baseVal})(-1);k.on("snap.util.getattr.fontSize",f)(-1);k.on("snap.util.getattr.font-size",f)(-1)});C.plugin(function(a,v,y,
M,A){function w(a){return a}function z(a){return function(b){return+b.toFixed(3)+a}}var d={"+":function(a,b){return a+b},"-":function(a,b){return a-b},"/":function(a,b){return a/b},"*":function(a,b){return a*b}},f=String,n=/[a-z]+$/i,u=/^\s*([+\-\/*])\s*=\s*([\d.eE+\-]+)\s*([^\d\s]+)?\s*$/;k.on("snap.util.attr",function(a){if(a=f(a).match(u)){var b=k.nt(),b=b.substring(b.lastIndexOf(".")+1),q=this.attr(b),e={};k.stop();var l=a[3]||"",r=q.match(n),s=d[a[1] ];r&&r==l?a=s(parseFloat(q),+a[2]):(q=this.asPX(b),
a=s(this.asPX(b),this.asPX(b,a[2]+l)));isNaN(q)||isNaN(a)||(e[b]=a,this.attr(e))}})(-10);k.on("snap.util.equal",function(a,b){var q=f(this.attr(a)||""),e=f(b).match(u);if(e){k.stop();var l=e[3]||"",r=q.match(n),s=d[e[1] ];if(r&&r==l)return{from:parseFloat(q),to:s(parseFloat(q),+e[2]),f:z(r)};q=this.asPX(a);return{from:q,to:s(q,this.asPX(a,e[2]+l)),f:w}}})(-10)});C.plugin(function(a,v,y,M,A){var w=y.prototype,z=a.is;w.rect=function(a,d,k,p,b,q){var e;null==q&&(q=b);z(a,"object")&&"[object Object]"==
a?e=a:null!=a&&(e={x:a,y:d,width:k,height:p},null!=b&&(e.rx=b,e.ry=q));return this.el("rect",e)};w.circle=function(a,d,k){var p;z(a,"object")&&"[object Object]"==a?p=a:null!=a&&(p={cx:a,cy:d,r:k});return this.el("circle",p)};var d=function(){function a(){this.parentNode.removeChild(this)}return function(d,k){var p=M.doc.createElement("img"),b=M.doc.body;p.style.cssText="position:absolute;left:-9999em;top:-9999em";p.onload=function(){k.call(p);p.onload=p.onerror=null;b.removeChild(p)};p.onerror=a;
b.appendChild(p);p.src=d}}();w.image=function(f,n,k,p,b){var q=this.el("image");if(z(f,"object")&&"src"in f)q.attr(f);else if(null!=f){var e={"xlink:href":f,preserveAspectRatio:"none"};null!=n&&null!=k&&(e.x=n,e.y=k);null!=p&&null!=b?(e.width=p,e.height=b):d(f,function(){a._.$(q.node,{width:this.offsetWidth,height:this.offsetHeight})});a._.$(q.node,e)}return q};w.ellipse=function(a,d,k,p){var b;z(a,"object")&&"[object Object]"==a?b=a:null!=a&&(b={cx:a,cy:d,rx:k,ry:p});return this.el("ellipse",b)};
w.path=function(a){var d;z(a,"object")&&!z(a,"array")?d=a:a&&(d={d:a});return this.el("path",d)};w.group=w.g=function(a){var d=this.el("g");1==arguments.length&&a&&!a.type?d.attr(a):arguments.length&&d.add(Array.prototype.slice.call(arguments,0));return d};w.svg=function(a,d,k,p,b,q,e,l){var r={};z(a,"object")&&null==d?r=a:(null!=a&&(r.x=a),null!=d&&(r.y=d),null!=k&&(r.width=k),null!=p&&(r.height=p),null!=b&&null!=q&&null!=e&&null!=l&&(r.viewBox=[b,q,e,l]));return this.el("svg",r)};w.mask=function(a){var d=
this.el("mask");1==arguments.length&&a&&!a.type?d.attr(a):arguments.length&&d.add(Array.prototype.slice.call(arguments,0));return d};w.ptrn=function(a,d,k,p,b,q,e,l){if(z(a,"object"))var r=a;else arguments.length?(r={},null!=a&&(r.x=a),null!=d&&(r.y=d),null!=k&&(r.width=k),null!=p&&(r.height=p),null!=b&&null!=q&&null!=e&&null!=l&&(r.viewBox=[b,q,e,l])):r={patternUnits:"userSpaceOnUse"};return this.el("pattern",r)};w.use=function(a){return null!=a?(make("use",this.node),a instanceof v&&(a.attr("id")||
a.attr({id:ID()}),a=a.attr("id")),this.el("use",{"xlink:href":a})):v.prototype.use.call(this)};w.text=function(a,d,k){var p={};z(a,"object")?p=a:null!=a&&(p={x:a,y:d,text:k||""});return this.el("text",p)};w.line=function(a,d,k,p){var b={};z(a,"object")?b=a:null!=a&&(b={x1:a,x2:k,y1:d,y2:p});return this.el("line",b)};w.polyline=function(a){1<arguments.length&&(a=Array.prototype.slice.call(arguments,0));var d={};z(a,"object")&&!z(a,"array")?d=a:null!=a&&(d={points:a});return this.el("polyline",d)};
w.polygon=function(a){1<arguments.length&&(a=Array.prototype.slice.call(arguments,0));var d={};z(a,"object")&&!z(a,"array")?d=a:null!=a&&(d={points:a});return this.el("polygon",d)};(function(){function d(){return this.selectAll("stop")}function n(b,d){var f=e("stop"),k={offset:+d+"%"};b=a.color(b);k["stop-color"]=b.hex;1>b.opacity&&(k["stop-opacity"]=b.opacity);e(f,k);this.node.appendChild(f);return this}function u(){if("linearGradient"==this.type){var b=e(this.node,"x1")||0,d=e(this.node,"x2")||
1,f=e(this.node,"y1")||0,k=e(this.node,"y2")||0;return a._.box(b,f,math.abs(d-b),math.abs(k-f))}b=this.node.r||0;return a._.box((this.node.cx||0.5)-b,(this.node.cy||0.5)-b,2*b,2*b)}function p(a,d){function f(a,b){for(var d=(b-u)/(a-w),e=w;e<a;e++)h[e].offset=+(+u+d*(e-w)).toFixed(2);w=a;u=b}var n=k("snap.util.grad.parse",null,d).firstDefined(),p;if(!n)return null;n.params.unshift(a);p="l"==n.type.toLowerCase()?b.apply(0,n.params):q.apply(0,n.params);n.type!=n.type.toLowerCase()&&e(p.node,{gradientUnits:"userSpaceOnUse"});
var h=n.stops,n=h.length,u=0,w=0;n--;for(var v=0;v<n;v++)"offset"in h[v]&&f(v,h[v].offset);h[n].offset=h[n].offset||100;f(n,h[n].offset);for(v=0;v<=n;v++){var y=h[v];p.addStop(y.color,y.offset)}return p}function b(b,k,p,q,w){b=a._.make("linearGradient",b);b.stops=d;b.addStop=n;b.getBBox=u;null!=k&&e(b.node,{x1:k,y1:p,x2:q,y2:w});return b}function q(b,k,p,q,w,h){b=a._.make("radialGradient",b);b.stops=d;b.addStop=n;b.getBBox=u;null!=k&&e(b.node,{cx:k,cy:p,r:q});null!=w&&null!=h&&e(b.node,{fx:w,fy:h});
return b}var e=a._.$;w.gradient=function(a){return p(this.defs,a)};w.gradientLinear=function(a,d,e,f){return b(this.defs,a,d,e,f)};w.gradientRadial=function(a,b,d,e,f){return q(this.defs,a,b,d,e,f)};w.toString=function(){var b=this.node.ownerDocument,d=b.createDocumentFragment(),b=b.createElement("div"),e=this.node.cloneNode(!0);d.appendChild(b);b.appendChild(e);a._.$(e,{xmlns:"http://www.w3.org/2000/svg"});b=b.innerHTML;d.removeChild(d.firstChild);return b};w.clear=function(){for(var a=this.node.firstChild,
b;a;)b=a.nextSibling,"defs"!=a.tagName?a.parentNode.removeChild(a):w.clear.call({node:a}),a=b}})()});C.plugin(function(a,k,y,M){function A(a){var b=A.ps=A.ps||{};b[a]?b[a].sleep=100:b[a]={sleep:100};setTimeout(function(){for(var d in b)b[L](d)&&d!=a&&(b[d].sleep--,!b[d].sleep&&delete b[d])});return b[a]}function w(a,b,d,e){null==a&&(a=b=d=e=0);null==b&&(b=a.y,d=a.width,e=a.height,a=a.x);return{x:a,y:b,width:d,w:d,height:e,h:e,x2:a+d,y2:b+e,cx:a+d/2,cy:b+e/2,r1:F.min(d,e)/2,r2:F.max(d,e)/2,r0:F.sqrt(d*
d+e*e)/2,path:s(a,b,d,e),vb:[a,b,d,e].join(" ")}}function z(){return this.join(",").replace(N,"$1")}function d(a){a=C(a);a.toString=z;return a}function f(a,b,d,h,f,k,l,n,p){if(null==p)return e(a,b,d,h,f,k,l,n);if(0>p||e(a,b,d,h,f,k,l,n)<p)p=void 0;else{var q=0.5,O=1-q,s;for(s=e(a,b,d,h,f,k,l,n,O);0.01<Z(s-p);)q/=2,O+=(s<p?1:-1)*q,s=e(a,b,d,h,f,k,l,n,O);p=O}return u(a,b,d,h,f,k,l,n,p)}function n(b,d){function e(a){return+(+a).toFixed(3)}return a._.cacher(function(a,h,l){a instanceof k&&(a=a.attr("d"));
a=I(a);for(var n,p,D,q,O="",s={},c=0,t=0,r=a.length;t<r;t++){D=a[t];if("M"==D[0])n=+D[1],p=+D[2];else{q=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6]);if(c+q>h){if(d&&!s.start){n=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6],h-c);O+=["C"+e(n.start.x),e(n.start.y),e(n.m.x),e(n.m.y),e(n.x),e(n.y)];if(l)return O;s.start=O;O=["M"+e(n.x),e(n.y)+"C"+e(n.n.x),e(n.n.y),e(n.end.x),e(n.end.y),e(D[5]),e(D[6])].join();c+=q;n=+D[5];p=+D[6];continue}if(!b&&!d)return n=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6],h-c)}c+=q;n=+D[5];p=+D[6]}O+=
D.shift()+D}s.end=O;return n=b?c:d?s:u(n,p,D[0],D[1],D[2],D[3],D[4],D[5],1)},null,a._.clone)}function u(a,b,d,e,h,f,k,l,n){var p=1-n,q=ma(p,3),s=ma(p,2),c=n*n,t=c*n,r=q*a+3*s*n*d+3*p*n*n*h+t*k,q=q*b+3*s*n*e+3*p*n*n*f+t*l,s=a+2*n*(d-a)+c*(h-2*d+a),t=b+2*n*(e-b)+c*(f-2*e+b),x=d+2*n*(h-d)+c*(k-2*h+d),c=e+2*n*(f-e)+c*(l-2*f+e);a=p*a+n*d;b=p*b+n*e;h=p*h+n*k;f=p*f+n*l;l=90-180*F.atan2(s-x,t-c)/S;return{x:r,y:q,m:{x:s,y:t},n:{x:x,y:c},start:{x:a,y:b},end:{x:h,y:f},alpha:l}}function p(b,d,e,h,f,n,k,l){a.is(b,
"array")||(b=[b,d,e,h,f,n,k,l]);b=U.apply(null,b);return w(b.min.x,b.min.y,b.max.x-b.min.x,b.max.y-b.min.y)}function b(a,b,d){return b>=a.x&&b<=a.x+a.width&&d>=a.y&&d<=a.y+a.height}function q(a,d){a=w(a);d=w(d);return b(d,a.x,a.y)||b(d,a.x2,a.y)||b(d,a.x,a.y2)||b(d,a.x2,a.y2)||b(a,d.x,d.y)||b(a,d.x2,d.y)||b(a,d.x,d.y2)||b(a,d.x2,d.y2)||(a.x<d.x2&&a.x>d.x||d.x<a.x2&&d.x>a.x)&&(a.y<d.y2&&a.y>d.y||d.y<a.y2&&d.y>a.y)}function e(a,b,d,e,h,f,n,k,l){null==l&&(l=1);l=(1<l?1:0>l?0:l)/2;for(var p=[-0.1252,
0.1252,-0.3678,0.3678,-0.5873,0.5873,-0.7699,0.7699,-0.9041,0.9041,-0.9816,0.9816],q=[0.2491,0.2491,0.2335,0.2335,0.2032,0.2032,0.1601,0.1601,0.1069,0.1069,0.0472,0.0472],s=0,c=0;12>c;c++)var t=l*p[c]+l,r=t*(t*(-3*a+9*d-9*h+3*n)+6*a-12*d+6*h)-3*a+3*d,t=t*(t*(-3*b+9*e-9*f+3*k)+6*b-12*e+6*f)-3*b+3*e,s=s+q[c]*F.sqrt(r*r+t*t);return l*s}function l(a,b,d){a=I(a);b=I(b);for(var h,f,l,n,k,s,r,O,x,c,t=d?0:[],w=0,v=a.length;w<v;w++)if(x=a[w],"M"==x[0])h=k=x[1],f=s=x[2];else{"C"==x[0]?(x=[h,f].concat(x.slice(1)),
h=x[6],f=x[7]):(x=[h,f,h,f,k,s,k,s],h=k,f=s);for(var G=0,y=b.length;G<y;G++)if(c=b[G],"M"==c[0])l=r=c[1],n=O=c[2];else{"C"==c[0]?(c=[l,n].concat(c.slice(1)),l=c[6],n=c[7]):(c=[l,n,l,n,r,O,r,O],l=r,n=O);var z;var K=x,B=c;z=d;var H=p(K),J=p(B);if(q(H,J)){for(var H=e.apply(0,K),J=e.apply(0,B),H=~~(H/8),J=~~(J/8),U=[],A=[],F={},M=z?0:[],P=0;P<H+1;P++){var C=u.apply(0,K.concat(P/H));U.push({x:C.x,y:C.y,t:P/H})}for(P=0;P<J+1;P++)C=u.apply(0,B.concat(P/J)),A.push({x:C.x,y:C.y,t:P/J});for(P=0;P<H;P++)for(K=
0;K<J;K++){var Q=U[P],L=U[P+1],B=A[K],C=A[K+1],N=0.001>Z(L.x-Q.x)?"y":"x",S=0.001>Z(C.x-B.x)?"y":"x",R;R=Q.x;var Y=Q.y,V=L.x,ea=L.y,fa=B.x,ga=B.y,ha=C.x,ia=C.y;if(W(R,V)<X(fa,ha)||X(R,V)>W(fa,ha)||W(Y,ea)<X(ga,ia)||X(Y,ea)>W(ga,ia))R=void 0;else{var $=(R*ea-Y*V)*(fa-ha)-(R-V)*(fa*ia-ga*ha),aa=(R*ea-Y*V)*(ga-ia)-(Y-ea)*(fa*ia-ga*ha),ja=(R-V)*(ga-ia)-(Y-ea)*(fa-ha);if(ja){var $=$/ja,aa=aa/ja,ja=+$.toFixed(2),ba=+aa.toFixed(2);R=ja<+X(R,V).toFixed(2)||ja>+W(R,V).toFixed(2)||ja<+X(fa,ha).toFixed(2)||
ja>+W(fa,ha).toFixed(2)||ba<+X(Y,ea).toFixed(2)||ba>+W(Y,ea).toFixed(2)||ba<+X(ga,ia).toFixed(2)||ba>+W(ga,ia).toFixed(2)?void 0:{x:$,y:aa}}else R=void 0}R&&F[R.x.toFixed(4)]!=R.y.toFixed(4)&&(F[R.x.toFixed(4)]=R.y.toFixed(4),Q=Q.t+Z((R[N]-Q[N])/(L[N]-Q[N]))*(L.t-Q.t),B=B.t+Z((R[S]-B[S])/(C[S]-B[S]))*(C.t-B.t),0<=Q&&1>=Q&&0<=B&&1>=B&&(z?M++:M.push({x:R.x,y:R.y,t1:Q,t2:B})))}z=M}else z=z?0:[];if(d)t+=z;else{H=0;for(J=z.length;H<J;H++)z[H].segment1=w,z[H].segment2=G,z[H].bez1=x,z[H].bez2=c;t=t.concat(z)}}}return t}
function r(a){var b=A(a);if(b.bbox)return C(b.bbox);if(!a)return w();a=I(a);for(var d=0,e=0,h=[],f=[],l,n=0,k=a.length;n<k;n++)l=a[n],"M"==l[0]?(d=l[1],e=l[2],h.push(d),f.push(e)):(d=U(d,e,l[1],l[2],l[3],l[4],l[5],l[6]),h=h.concat(d.min.x,d.max.x),f=f.concat(d.min.y,d.max.y),d=l[5],e=l[6]);a=X.apply(0,h);l=X.apply(0,f);h=W.apply(0,h);f=W.apply(0,f);f=w(a,l,h-a,f-l);b.bbox=C(f);return f}function s(a,b,d,e,h){if(h)return[["M",+a+ +h,b],["l",d-2*h,0],["a",h,h,0,0,1,h,h],["l",0,e-2*h],["a",h,h,0,0,1,
-h,h],["l",2*h-d,0],["a",h,h,0,0,1,-h,-h],["l",0,2*h-e],["a",h,h,0,0,1,h,-h],["z"] ];a=[["M",a,b],["l",d,0],["l",0,e],["l",-d,0],["z"] ];a.toString=z;return a}function x(a,b,d,e,h){null==h&&null==e&&(e=d);a=+a;b=+b;d=+d;e=+e;if(null!=h){var f=Math.PI/180,l=a+d*Math.cos(-e*f);a+=d*Math.cos(-h*f);var n=b+d*Math.sin(-e*f);b+=d*Math.sin(-h*f);d=[["M",l,n],["A",d,d,0,+(180<h-e),0,a,b] ]}else d=[["M",a,b],["m",0,-e],["a",d,e,0,1,1,0,2*e],["a",d,e,0,1,1,0,-2*e],["z"] ];d.toString=z;return d}function G(b){var e=
A(b);if(e.abs)return d(e.abs);Q(b,"array")&&Q(b&&b[0],"array")||(b=a.parsePathString(b));if(!b||!b.length)return[["M",0,0] ];var h=[],f=0,l=0,n=0,k=0,p=0;"M"==b[0][0]&&(f=+b[0][1],l=+b[0][2],n=f,k=l,p++,h[0]=["M",f,l]);for(var q=3==b.length&&"M"==b[0][0]&&"R"==b[1][0].toUpperCase()&&"Z"==b[2][0].toUpperCase(),s,r,w=p,c=b.length;w<c;w++){h.push(s=[]);r=b[w];p=r[0];if(p!=p.toUpperCase())switch(s[0]=p.toUpperCase(),s[0]){case "A":s[1]=r[1];s[2]=r[2];s[3]=r[3];s[4]=r[4];s[5]=r[5];s[6]=+r[6]+f;s[7]=+r[7]+
l;break;case "V":s[1]=+r[1]+l;break;case "H":s[1]=+r[1]+f;break;case "R":for(var t=[f,l].concat(r.slice(1)),u=2,v=t.length;u<v;u++)t[u]=+t[u]+f,t[++u]=+t[u]+l;h.pop();h=h.concat(P(t,q));break;case "O":h.pop();t=x(f,l,r[1],r[2]);t.push(t[0]);h=h.concat(t);break;case "U":h.pop();h=h.concat(x(f,l,r[1],r[2],r[3]));s=["U"].concat(h[h.length-1].slice(-2));break;case "M":n=+r[1]+f,k=+r[2]+l;default:for(u=1,v=r.length;u<v;u++)s[u]=+r[u]+(u%2?f:l)}else if("R"==p)t=[f,l].concat(r.slice(1)),h.pop(),h=h.concat(P(t,
q)),s=["R"].concat(r.slice(-2));else if("O"==p)h.pop(),t=x(f,l,r[1],r[2]),t.push(t[0]),h=h.concat(t);else if("U"==p)h.pop(),h=h.concat(x(f,l,r[1],r[2],r[3])),s=["U"].concat(h[h.length-1].slice(-2));else for(t=0,u=r.length;t<u;t++)s[t]=r[t];p=p.toUpperCase();if("O"!=p)switch(s[0]){case "Z":f=+n;l=+k;break;case "H":f=s[1];break;case "V":l=s[1];break;case "M":n=s[s.length-2],k=s[s.length-1];default:f=s[s.length-2],l=s[s.length-1]}}h.toString=z;e.abs=d(h);return h}function h(a,b,d,e){return[a,b,d,e,d,
e]}function J(a,b,d,e,h,f){var l=1/3,n=2/3;return[l*a+n*d,l*b+n*e,l*h+n*d,l*f+n*e,h,f]}function K(b,d,e,h,f,l,n,k,p,s){var r=120*S/180,q=S/180*(+f||0),c=[],t,x=a._.cacher(function(a,b,c){var d=a*F.cos(c)-b*F.sin(c);a=a*F.sin(c)+b*F.cos(c);return{x:d,y:a}});if(s)v=s[0],t=s[1],l=s[2],u=s[3];else{t=x(b,d,-q);b=t.x;d=t.y;t=x(k,p,-q);k=t.x;p=t.y;F.cos(S/180*f);F.sin(S/180*f);t=(b-k)/2;v=(d-p)/2;u=t*t/(e*e)+v*v/(h*h);1<u&&(u=F.sqrt(u),e*=u,h*=u);var u=e*e,w=h*h,u=(l==n?-1:1)*F.sqrt(Z((u*w-u*v*v-w*t*t)/
(u*v*v+w*t*t)));l=u*e*v/h+(b+k)/2;var u=u*-h*t/e+(d+p)/2,v=F.asin(((d-u)/h).toFixed(9));t=F.asin(((p-u)/h).toFixed(9));v=b<l?S-v:v;t=k<l?S-t:t;0>v&&(v=2*S+v);0>t&&(t=2*S+t);n&&v>t&&(v-=2*S);!n&&t>v&&(t-=2*S)}if(Z(t-v)>r){var c=t,w=k,G=p;t=v+r*(n&&t>v?1:-1);k=l+e*F.cos(t);p=u+h*F.sin(t);c=K(k,p,e,h,f,0,n,w,G,[t,c,l,u])}l=t-v;f=F.cos(v);r=F.sin(v);n=F.cos(t);t=F.sin(t);l=F.tan(l/4);e=4/3*e*l;l*=4/3*h;h=[b,d];b=[b+e*r,d-l*f];d=[k+e*t,p-l*n];k=[k,p];b[0]=2*h[0]-b[0];b[1]=2*h[1]-b[1];if(s)return[b,d,k].concat(c);
c=[b,d,k].concat(c).join().split(",");s=[];k=0;for(p=c.length;k<p;k++)s[k]=k%2?x(c[k-1],c[k],q).y:x(c[k],c[k+1],q).x;return s}function U(a,b,d,e,h,f,l,k){for(var n=[],p=[[],[] ],s,r,c,t,q=0;2>q;++q)0==q?(r=6*a-12*d+6*h,s=-3*a+9*d-9*h+3*l,c=3*d-3*a):(r=6*b-12*e+6*f,s=-3*b+9*e-9*f+3*k,c=3*e-3*b),1E-12>Z(s)?1E-12>Z(r)||(s=-c/r,0<s&&1>s&&n.push(s)):(t=r*r-4*c*s,c=F.sqrt(t),0>t||(t=(-r+c)/(2*s),0<t&&1>t&&n.push(t),s=(-r-c)/(2*s),0<s&&1>s&&n.push(s)));for(r=q=n.length;q--;)s=n[q],c=1-s,p[0][q]=c*c*c*a+3*
c*c*s*d+3*c*s*s*h+s*s*s*l,p[1][q]=c*c*c*b+3*c*c*s*e+3*c*s*s*f+s*s*s*k;p[0][r]=a;p[1][r]=b;p[0][r+1]=l;p[1][r+1]=k;p[0].length=p[1].length=r+2;return{min:{x:X.apply(0,p[0]),y:X.apply(0,p[1])},max:{x:W.apply(0,p[0]),y:W.apply(0,p[1])}}}function I(a,b){var e=!b&&A(a);if(!b&&e.curve)return d(e.curve);var f=G(a),l=b&&G(b),n={x:0,y:0,bx:0,by:0,X:0,Y:0,qx:null,qy:null},k={x:0,y:0,bx:0,by:0,X:0,Y:0,qx:null,qy:null},p=function(a,b,c){if(!a)return["C",b.x,b.y,b.x,b.y,b.x,b.y];a[0]in{T:1,Q:1}||(b.qx=b.qy=null);
switch(a[0]){case "M":b.X=a[1];b.Y=a[2];break;case "A":a=["C"].concat(K.apply(0,[b.x,b.y].concat(a.slice(1))));break;case "S":"C"==c||"S"==c?(c=2*b.x-b.bx,b=2*b.y-b.by):(c=b.x,b=b.y);a=["C",c,b].concat(a.slice(1));break;case "T":"Q"==c||"T"==c?(b.qx=2*b.x-b.qx,b.qy=2*b.y-b.qy):(b.qx=b.x,b.qy=b.y);a=["C"].concat(J(b.x,b.y,b.qx,b.qy,a[1],a[2]));break;case "Q":b.qx=a[1];b.qy=a[2];a=["C"].concat(J(b.x,b.y,a[1],a[2],a[3],a[4]));break;case "L":a=["C"].concat(h(b.x,b.y,a[1],a[2]));break;case "H":a=["C"].concat(h(b.x,
b.y,a[1],b.y));break;case "V":a=["C"].concat(h(b.x,b.y,b.x,a[1]));break;case "Z":a=["C"].concat(h(b.x,b.y,b.X,b.Y))}return a},s=function(a,b){if(7<a[b].length){a[b].shift();for(var c=a[b];c.length;)q[b]="A",l&&(u[b]="A"),a.splice(b++,0,["C"].concat(c.splice(0,6)));a.splice(b,1);v=W(f.length,l&&l.length||0)}},r=function(a,b,c,d,e){a&&b&&"M"==a[e][0]&&"M"!=b[e][0]&&(b.splice(e,0,["M",d.x,d.y]),c.bx=0,c.by=0,c.x=a[e][1],c.y=a[e][2],v=W(f.length,l&&l.length||0))},q=[],u=[],c="",t="",x=0,v=W(f.length,
l&&l.length||0);for(;x<v;x++){f[x]&&(c=f[x][0]);"C"!=c&&(q[x]=c,x&&(t=q[x-1]));f[x]=p(f[x],n,t);"A"!=q[x]&&"C"==c&&(q[x]="C");s(f,x);l&&(l[x]&&(c=l[x][0]),"C"!=c&&(u[x]=c,x&&(t=u[x-1])),l[x]=p(l[x],k,t),"A"!=u[x]&&"C"==c&&(u[x]="C"),s(l,x));r(f,l,n,k,x);r(l,f,k,n,x);var w=f[x],z=l&&l[x],y=w.length,U=l&&z.length;n.x=w[y-2];n.y=w[y-1];n.bx=$(w[y-4])||n.x;n.by=$(w[y-3])||n.y;k.bx=l&&($(z[U-4])||k.x);k.by=l&&($(z[U-3])||k.y);k.x=l&&z[U-2];k.y=l&&z[U-1]}l||(e.curve=d(f));return l?[f,l]:f}function P(a,
b){for(var d=[],e=0,h=a.length;h-2*!b>e;e+=2){var f=[{x:+a[e-2],y:+a[e-1]},{x:+a[e],y:+a[e+1]},{x:+a[e+2],y:+a[e+3]},{x:+a[e+4],y:+a[e+5]}];b?e?h-4==e?f[3]={x:+a[0],y:+a[1]}:h-2==e&&(f[2]={x:+a[0],y:+a[1]},f[3]={x:+a[2],y:+a[3]}):f[0]={x:+a[h-2],y:+a[h-1]}:h-4==e?f[3]=f[2]:e||(f[0]={x:+a[e],y:+a[e+1]});d.push(["C",(-f[0].x+6*f[1].x+f[2].x)/6,(-f[0].y+6*f[1].y+f[2].y)/6,(f[1].x+6*f[2].x-f[3].x)/6,(f[1].y+6*f[2].y-f[3].y)/6,f[2].x,f[2].y])}return d}y=k.prototype;var Q=a.is,C=a._.clone,L="hasOwnProperty",
N=/,?([a-z]),?/gi,$=parseFloat,F=Math,S=F.PI,X=F.min,W=F.max,ma=F.pow,Z=F.abs;M=n(1);var na=n(),ba=n(0,1),V=a._unit2px;a.path=A;a.path.getTotalLength=M;a.path.getPointAtLength=na;a.path.getSubpath=function(a,b,d){if(1E-6>this.getTotalLength(a)-d)return ba(a,b).end;a=ba(a,d,1);return b?ba(a,b).end:a};y.getTotalLength=function(){if(this.node.getTotalLength)return this.node.getTotalLength()};y.getPointAtLength=function(a){return na(this.attr("d"),a)};y.getSubpath=function(b,d){return a.path.getSubpath(this.attr("d"),
b,d)};a._.box=w;a.path.findDotsAtSegment=u;a.path.bezierBBox=p;a.path.isPointInsideBBox=b;a.path.isBBoxIntersect=q;a.path.intersection=function(a,b){return l(a,b)};a.path.intersectionNumber=function(a,b){return l(a,b,1)};a.path.isPointInside=function(a,d,e){var h=r(a);return b(h,d,e)&&1==l(a,[["M",d,e],["H",h.x2+10] ],1)%2};a.path.getBBox=r;a.path.get={path:function(a){return a.attr("path")},circle:function(a){a=V(a);return x(a.cx,a.cy,a.r)},ellipse:function(a){a=V(a);return x(a.cx||0,a.cy||0,a.rx,
a.ry)},rect:function(a){a=V(a);return s(a.x||0,a.y||0,a.width,a.height,a.rx,a.ry)},image:function(a){a=V(a);return s(a.x||0,a.y||0,a.width,a.height)},line:function(a){return"M"+[a.attr("x1")||0,a.attr("y1")||0,a.attr("x2"),a.attr("y2")]},polyline:function(a){return"M"+a.attr("points")},polygon:function(a){return"M"+a.attr("points")+"z"},deflt:function(a){a=a.node.getBBox();return s(a.x,a.y,a.width,a.height)}};a.path.toRelative=function(b){var e=A(b),h=String.prototype.toLowerCase;if(e.rel)return d(e.rel);
a.is(b,"array")&&a.is(b&&b[0],"array")||(b=a.parsePathString(b));var f=[],l=0,n=0,k=0,p=0,s=0;"M"==b[0][0]&&(l=b[0][1],n=b[0][2],k=l,p=n,s++,f.push(["M",l,n]));for(var r=b.length;s<r;s++){var q=f[s]=[],x=b[s];if(x[0]!=h.call(x[0]))switch(q[0]=h.call(x[0]),q[0]){case "a":q[1]=x[1];q[2]=x[2];q[3]=x[3];q[4]=x[4];q[5]=x[5];q[6]=+(x[6]-l).toFixed(3);q[7]=+(x[7]-n).toFixed(3);break;case "v":q[1]=+(x[1]-n).toFixed(3);break;case "m":k=x[1],p=x[2];default:for(var c=1,t=x.length;c<t;c++)q[c]=+(x[c]-(c%2?l:
n)).toFixed(3)}else for(f[s]=[],"m"==x[0]&&(k=x[1]+l,p=x[2]+n),q=0,c=x.length;q<c;q++)f[s][q]=x[q];x=f[s].length;switch(f[s][0]){case "z":l=k;n=p;break;case "h":l+=+f[s][x-1];break;case "v":n+=+f[s][x-1];break;default:l+=+f[s][x-2],n+=+f[s][x-1]}}f.toString=z;e.rel=d(f);return f};a.path.toAbsolute=G;a.path.toCubic=I;a.path.map=function(a,b){if(!b)return a;var d,e,h,f,l,n,k;a=I(a);h=0;for(l=a.length;h<l;h++)for(k=a[h],f=1,n=k.length;f<n;f+=2)d=b.x(k[f],k[f+1]),e=b.y(k[f],k[f+1]),k[f]=d,k[f+1]=e;return a};
a.path.toString=z;a.path.clone=d});C.plugin(function(a,v,y,C){var A=Math.max,w=Math.min,z=function(a){this.items=[];this.bindings={};this.length=0;this.type="set";if(a)for(var f=0,n=a.length;f<n;f++)a[f]&&(this[this.items.length]=this.items[this.items.length]=a[f],this.length++)};v=z.prototype;v.push=function(){for(var a,f,n=0,k=arguments.length;n<k;n++)if(a=arguments[n])f=this.items.length,this[f]=this.items[f]=a,this.length++;return this};v.pop=function(){this.length&&delete this[this.length--];
return this.items.pop()};v.forEach=function(a,f){for(var n=0,k=this.items.length;n<k&&!1!==a.call(f,this.items[n],n);n++);return this};v.animate=function(d,f,n,u){"function"!=typeof n||n.length||(u=n,n=L.linear);d instanceof a._.Animation&&(u=d.callback,n=d.easing,f=n.dur,d=d.attr);var p=arguments;if(a.is(d,"array")&&a.is(p[p.length-1],"array"))var b=!0;var q,e=function(){q?this.b=q:q=this.b},l=0,r=u&&function(){l++==this.length&&u.call(this)};return this.forEach(function(a,l){k.once("snap.animcreated."+
a.id,e);b?p[l]&&a.animate.apply(a,p[l]):a.animate(d,f,n,r)})};v.remove=function(){for(;this.length;)this.pop().remove();return this};v.bind=function(a,f,k){var u={};if("function"==typeof f)this.bindings[a]=f;else{var p=k||a;this.bindings[a]=function(a){u[p]=a;f.attr(u)}}return this};v.attr=function(a){var f={},k;for(k in a)if(this.bindings[k])this.bindings[k](a[k]);else f[k]=a[k];a=0;for(k=this.items.length;a<k;a++)this.items[a].attr(f);return this};v.clear=function(){for(;this.length;)this.pop()};
v.splice=function(a,f,k){a=0>a?A(this.length+a,0):a;f=A(0,w(this.length-a,f));var u=[],p=[],b=[],q;for(q=2;q<arguments.length;q++)b.push(arguments[q]);for(q=0;q<f;q++)p.push(this[a+q]);for(;q<this.length-a;q++)u.push(this[a+q]);var e=b.length;for(q=0;q<e+u.length;q++)this.items[a+q]=this[a+q]=q<e?b[q]:u[q-e];for(q=this.items.length=this.length-=f-e;this[q];)delete this[q++];return new z(p)};v.exclude=function(a){for(var f=0,k=this.length;f<k;f++)if(this[f]==a)return this.splice(f,1),!0;return!1};
v.insertAfter=function(a){for(var f=this.items.length;f--;)this.items[f].insertAfter(a);return this};v.getBBox=function(){for(var a=[],f=[],k=[],u=[],p=this.items.length;p--;)if(!this.items[p].removed){var b=this.items[p].getBBox();a.push(b.x);f.push(b.y);k.push(b.x+b.width);u.push(b.y+b.height)}a=w.apply(0,a);f=w.apply(0,f);k=A.apply(0,k);u=A.apply(0,u);return{x:a,y:f,x2:k,y2:u,width:k-a,height:u-f,cx:a+(k-a)/2,cy:f+(u-f)/2}};v.clone=function(a){a=new z;for(var f=0,k=this.items.length;f<k;f++)a.push(this.items[f].clone());
return a};v.toString=function(){return"Snap\u2018s set"};v.type="set";a.set=function(){var a=new z;arguments.length&&a.push.apply(a,Array.prototype.slice.call(arguments,0));return a}});C.plugin(function(a,v,y,C){function A(a){var b=a[0];switch(b.toLowerCase()){case "t":return[b,0,0];case "m":return[b,1,0,0,1,0,0];case "r":return 4==a.length?[b,0,a[2],a[3] ]:[b,0];case "s":return 5==a.length?[b,1,1,a[3],a[4] ]:3==a.length?[b,1,1]:[b,1]}}function w(b,d,f){d=q(d).replace(/\.{3}|\u2026/g,b);b=a.parseTransformString(b)||
[];d=a.parseTransformString(d)||[];for(var k=Math.max(b.length,d.length),p=[],v=[],h=0,w,z,y,I;h<k;h++){y=b[h]||A(d[h]);I=d[h]||A(y);if(y[0]!=I[0]||"r"==y[0].toLowerCase()&&(y[2]!=I[2]||y[3]!=I[3])||"s"==y[0].toLowerCase()&&(y[3]!=I[3]||y[4]!=I[4])){b=a._.transform2matrix(b,f());d=a._.transform2matrix(d,f());p=[["m",b.a,b.b,b.c,b.d,b.e,b.f] ];v=[["m",d.a,d.b,d.c,d.d,d.e,d.f] ];break}p[h]=[];v[h]=[];w=0;for(z=Math.max(y.length,I.length);w<z;w++)w in y&&(p[h][w]=y[w]),w in I&&(v[h][w]=I[w])}return{from:u(p),
to:u(v),f:n(p)}}function z(a){return a}function d(a){return function(b){return+b.toFixed(3)+a}}function f(b){return a.rgb(b[0],b[1],b[2])}function n(a){var b=0,d,f,k,n,h,p,q=[];d=0;for(f=a.length;d<f;d++){h="[";p=['"'+a[d][0]+'"'];k=1;for(n=a[d].length;k<n;k++)p[k]="val["+b++ +"]";h+=p+"]";q[d]=h}return Function("val","return Snap.path.toString.call(["+q+"])")}function u(a){for(var b=[],d=0,f=a.length;d<f;d++)for(var k=1,n=a[d].length;k<n;k++)b.push(a[d][k]);return b}var p={},b=/[a-z]+$/i,q=String;
p.stroke=p.fill="colour";v.prototype.equal=function(a,b){return k("snap.util.equal",this,a,b).firstDefined()};k.on("snap.util.equal",function(e,k){var r,s;r=q(this.attr(e)||"");var x=this;if(r==+r&&k==+k)return{from:+r,to:+k,f:z};if("colour"==p[e])return r=a.color(r),s=a.color(k),{from:[r.r,r.g,r.b,r.opacity],to:[s.r,s.g,s.b,s.opacity],f:f};if("transform"==e||"gradientTransform"==e||"patternTransform"==e)return k instanceof a.Matrix&&(k=k.toTransformString()),a._.rgTransform.test(k)||(k=a._.svgTransform2string(k)),
w(r,k,function(){return x.getBBox(1)});if("d"==e||"path"==e)return r=a.path.toCubic(r,k),{from:u(r[0]),to:u(r[1]),f:n(r[0])};if("points"==e)return r=q(r).split(a._.separator),s=q(k).split(a._.separator),{from:r,to:s,f:function(a){return a}};aUnit=r.match(b);s=q(k).match(b);return aUnit&&aUnit==s?{from:parseFloat(r),to:parseFloat(k),f:d(aUnit)}:{from:this.asPX(e),to:this.asPX(e,k),f:z}})});C.plugin(function(a,v,y,C){var A=v.prototype,w="createTouch"in C.doc;v="click dblclick mousedown mousemove mouseout mouseover mouseup touchstart touchmove touchend touchcancel".split(" ");
var z={mousedown:"touchstart",mousemove:"touchmove",mouseup:"touchend"},d=function(a,b){var d="y"==a?"scrollTop":"scrollLeft",e=b&&b.node?b.node.ownerDocument:C.doc;return e[d in e.documentElement?"documentElement":"body"][d]},f=function(){this.returnValue=!1},n=function(){return this.originalEvent.preventDefault()},u=function(){this.cancelBubble=!0},p=function(){return this.originalEvent.stopPropagation()},b=function(){if(C.doc.addEventListener)return function(a,b,e,f){var k=w&&z[b]?z[b]:b,l=function(k){var l=
d("y",f),q=d("x",f);if(w&&z.hasOwnProperty(b))for(var r=0,u=k.targetTouches&&k.targetTouches.length;r<u;r++)if(k.targetTouches[r].target==a||a.contains(k.targetTouches[r].target)){u=k;k=k.targetTouches[r];k.originalEvent=u;k.preventDefault=n;k.stopPropagation=p;break}return e.call(f,k,k.clientX+q,k.clientY+l)};b!==k&&a.addEventListener(b,l,!1);a.addEventListener(k,l,!1);return function(){b!==k&&a.removeEventListener(b,l,!1);a.removeEventListener(k,l,!1);return!0}};if(C.doc.attachEvent)return function(a,
b,e,h){var k=function(a){a=a||h.node.ownerDocument.window.event;var b=d("y",h),k=d("x",h),k=a.clientX+k,b=a.clientY+b;a.preventDefault=a.preventDefault||f;a.stopPropagation=a.stopPropagation||u;return e.call(h,a,k,b)};a.attachEvent("on"+b,k);return function(){a.detachEvent("on"+b,k);return!0}}}(),q=[],e=function(a){for(var b=a.clientX,e=a.clientY,f=d("y"),l=d("x"),n,p=q.length;p--;){n=q[p];if(w)for(var r=a.touches&&a.touches.length,u;r--;){if(u=a.touches[r],u.identifier==n.el._drag.id||n.el.node.contains(u.target)){b=
u.clientX;e=u.clientY;(a.originalEvent?a.originalEvent:a).preventDefault();break}}else a.preventDefault();b+=l;e+=f;k("snap.drag.move."+n.el.id,n.move_scope||n.el,b-n.el._drag.x,e-n.el._drag.y,b,e,a)}},l=function(b){a.unmousemove(e).unmouseup(l);for(var d=q.length,f;d--;)f=q[d],f.el._drag={},k("snap.drag.end."+f.el.id,f.end_scope||f.start_scope||f.move_scope||f.el,b);q=[]};for(y=v.length;y--;)(function(d){a[d]=A[d]=function(e,f){a.is(e,"function")&&(this.events=this.events||[],this.events.push({name:d,
f:e,unbind:b(this.node||document,d,e,f||this)}));return this};a["un"+d]=A["un"+d]=function(a){for(var b=this.events||[],e=b.length;e--;)if(b[e].name==d&&(b[e].f==a||!a)){b[e].unbind();b.splice(e,1);!b.length&&delete this.events;break}return this}})(v[y]);A.hover=function(a,b,d,e){return this.mouseover(a,d).mouseout(b,e||d)};A.unhover=function(a,b){return this.unmouseover(a).unmouseout(b)};var r=[];A.drag=function(b,d,f,h,n,p){function u(r,v,w){(r.originalEvent||r).preventDefault();this._drag.x=v;
this._drag.y=w;this._drag.id=r.identifier;!q.length&&a.mousemove(e).mouseup(l);q.push({el:this,move_scope:h,start_scope:n,end_scope:p});d&&k.on("snap.drag.start."+this.id,d);b&&k.on("snap.drag.move."+this.id,b);f&&k.on("snap.drag.end."+this.id,f);k("snap.drag.start."+this.id,n||h||this,v,w,r)}if(!arguments.length){var v;return this.drag(function(a,b){this.attr({transform:v+(v?"T":"t")+[a,b]})},function(){v=this.transform().local})}this._drag={};r.push({el:this,start:u});this.mousedown(u);return this};
A.undrag=function(){for(var b=r.length;b--;)r[b].el==this&&(this.unmousedown(r[b].start),r.splice(b,1),k.unbind("snap.drag.*."+this.id));!r.length&&a.unmousemove(e).unmouseup(l);return this}});C.plugin(function(a,v,y,C){y=y.prototype;var A=/^\s*url\((.+)\)/,w=String,z=a._.$;a.filter={};y.filter=function(d){var f=this;"svg"!=f.type&&(f=f.paper);d=a.parse(w(d));var k=a._.id(),u=z("filter");z(u,{id:k,filterUnits:"userSpaceOnUse"});u.appendChild(d.node);f.defs.appendChild(u);return new v(u)};k.on("snap.util.getattr.filter",
function(){k.stop();var d=z(this.node,"filter");if(d)return(d=w(d).match(A))&&a.select(d[1])});k.on("snap.util.attr.filter",function(d){if(d instanceof v&&"filter"==d.type){k.stop();var f=d.node.id;f||(z(d.node,{id:d.id}),f=d.id);z(this.node,{filter:a.url(f)})}d&&"none"!=d||(k.stop(),this.node.removeAttribute("filter"))});a.filter.blur=function(d,f){null==d&&(d=2);return a.format('<feGaussianBlur stdDeviation="{def}"/>',{def:null==f?d:[d,f]})};a.filter.blur.toString=function(){return this()};a.filter.shadow=
function(d,f,k,u,p){"string"==typeof k&&(p=u=k,k=4);"string"!=typeof u&&(p=u,u="#000");null==k&&(k=4);null==p&&(p=1);null==d&&(d=0,f=2);null==f&&(f=d);u=a.color(u||"#000");return a.format('<feGaussianBlur in="SourceAlpha" stdDeviation="{blur}"/><feOffset dx="{dx}" dy="{dy}" result="offsetblur"/><feFlood flood-color="{color}"/><feComposite in2="offsetblur" operator="in"/><feComponentTransfer><feFuncA type="linear" slope="{opacity}"/></feComponentTransfer><feMerge><feMergeNode/><feMergeNode in="SourceGraphic"/></feMerge>',
{color:u,dx:d,dy:f,blur:k,opacity:p})};a.filter.shadow.toString=function(){return this()};a.filter.grayscale=function(d){null==d&&(d=1);return a.format('<feColorMatrix type="matrix" values="{a} {b} {c} 0 0 {d} {e} {f} 0 0 {g} {b} {h} 0 0 0 0 0 1 0"/>',{a:0.2126+0.7874*(1-d),b:0.7152-0.7152*(1-d),c:0.0722-0.0722*(1-d),d:0.2126-0.2126*(1-d),e:0.7152+0.2848*(1-d),f:0.0722-0.0722*(1-d),g:0.2126-0.2126*(1-d),h:0.0722+0.9278*(1-d)})};a.filter.grayscale.toString=function(){return this()};a.filter.sepia=
function(d){null==d&&(d=1);return a.format('<feColorMatrix type="matrix" values="{a} {b} {c} 0 0 {d} {e} {f} 0 0 {g} {h} {i} 0 0 0 0 0 1 0"/>',{a:0.393+0.607*(1-d),b:0.769-0.769*(1-d),c:0.189-0.189*(1-d),d:0.349-0.349*(1-d),e:0.686+0.314*(1-d),f:0.168-0.168*(1-d),g:0.272-0.272*(1-d),h:0.534-0.534*(1-d),i:0.131+0.869*(1-d)})};a.filter.sepia.toString=function(){return this()};a.filter.saturate=function(d){null==d&&(d=1);return a.format('<feColorMatrix type="saturate" values="{amount}"/>',{amount:1-
d})};a.filter.saturate.toString=function(){return this()};a.filter.hueRotate=function(d){return a.format('<feColorMatrix type="hueRotate" values="{angle}"/>',{angle:d||0})};a.filter.hueRotate.toString=function(){return this()};a.filter.invert=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="table" tableValues="{amount} {amount2}"/><feFuncG type="table" tableValues="{amount} {amount2}"/><feFuncB type="table" tableValues="{amount} {amount2}"/></feComponentTransfer>',{amount:d,
amount2:1-d})};a.filter.invert.toString=function(){return this()};a.filter.brightness=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="linear" slope="{amount}"/><feFuncG type="linear" slope="{amount}"/><feFuncB type="linear" slope="{amount}"/></feComponentTransfer>',{amount:d})};a.filter.brightness.toString=function(){return this()};a.filter.contrast=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="linear" slope="{amount}" intercept="{amount2}"/><feFuncG type="linear" slope="{amount}" intercept="{amount2}"/><feFuncB type="linear" slope="{amount}" intercept="{amount2}"/></feComponentTransfer>',
{amount:d,amount2:0.5-d/2})};a.filter.contrast.toString=function(){return this()}});return C});

]]> </script>
</svg>



<div class="markdown"><p>We can verify the correctness by checking the operator fidelity</p>
</div>

<pre class='language-julia'><code class='language-julia'>operator_fidelity(qft3' * qft3, igate(3))</code></pre>
<pre id='var-hash859977' class='code-output documenter-example-output'>0.9999999999999994</pre>


<div class="markdown"><p>Operator fidelity is defined as</p>
<p class="tex">$$F&#40;A, B&#41; &#61; \frac&#123;1&#125;&#123;D&#125;|&#123;\rm Tr&#125;&#40;A, B&#41;|$$</p>
<p>where <span class="tex">$D$</span> is the Hiltert space dimension.</p>
</div>


<div class="markdown"><ol start="2">
<li><p>Reflexive: <span class="tex">$A^2 &#61; \mathbb&#123;1&#125;$</span></p>
</li>
</ol>
</div>

<pre class='language-julia'><code class='language-julia'>isreflexive(X)</code></pre>
<pre id='var-hash878624' class='code-output documenter-example-output'>true</pre>


<div class="markdown"><p class="tex">$$&#123;\rm rot&#125;&#40;G, \theta&#41; &#61; \cos \frac&#123;\theta&#125;&#123;2&#125;I - i \sin \frac&#123;\theta&#125;&#123;2&#125; G$$</p>
</div>


<div class="markdown"><p><span class="tex">$\exp&#40;iG\theta/2&#41; &#61; \uparrow$</span></p>
</div>

<pre class='language-julia'><code class='language-julia'>mat(rot(X, 0.4))</code></pre>
<pre id='var-hash300348' class='code-output documenter-example-output'>2×2 Matrix{ComplexF64}:
 0.980067+0.0im            0.0-0.198669im
      0.0-0.198669im  0.980067+0.0im</pre>


<div class="markdown"><p>The Heisenberg Hamiltonian is hermitian</p>
</div>


<div class="markdown"><p>Only reflexive operators are allowed in rotation gate as the generator. The following code will error:</p>
</div>

<pre class='language-julia'><code class='language-julia'># rot(X + Y, 0.5)</code></pre>


<pre class='language-julia'><code class='language-julia'>X + Y</code></pre>
<pre id='var-hash103185' class='code-output documenter-example-output'>nqubits: 1
+
├─ X
└─ Y</pre>


<div class="markdown"><ol start="3">
<li><p>Unitary: <span class="tex">$A A^\dagger &#61; \mathbb&#123;1&#125;$</span></p>
</li>
</ol>
</div>

<pre class='language-julia'><code class='language-julia'>isunitary(X)</code></pre>
<pre id='var-hash129972' class='code-output documenter-example-output'>true</pre>


<div class="markdown"><p>Unitary is required by the back propagation algorithm &#40;will not be covered today&#41;.</p>
</div>


<div class="markdown"><ol start="4">
<li><p>Commutative: <span class="tex">$&#91;A, B&#93; &#61; \mathbb&#123;0&#125;$</span></p>
</li>
</ol>
</div>

<pre class='language-julia'><code class='language-julia'>iscommute(X, Y)</code></pre>
<pre id='var-hash172120' class='code-output documenter-example-output'>false</pre>


<div class="markdown"><p>These properties does not require calling the <code>mat</code>, unless it falls back.</p>
</div>

<pre class='language-julia'><code class='language-julia'>iscommute(put(1000, 4=&gt;X), put(1000, 6=&gt;Y))</code></pre>
<pre id='var-hash507660' class='code-output documenter-example-output'>true</pre>

<pre class='language-julia'><code class='language-julia'>iscommute(put(1000, 4=&gt;X), put(1000, 4=&gt;Y))</code></pre>
<pre id='var-hash517530' class='code-output documenter-example-output'>false</pre>

<pre class='language-julia'><code class='language-julia'>ishermitian(EasyBuild.heisenberg(1000))</code></pre>
<pre id='var-hash109640' class='code-output documenter-example-output'>true</pre>

<pre class='language-julia'><code class='language-julia'>EasyBuild.heisenberg(3)</code></pre>
<pre id='var-hash180806' class='code-output documenter-example-output'>nqubits: 3
+
├─ +
│  ├─ repeat on (1, 2)
│  │  └─ X
│  ├─ repeat on (1, 2)
│  │  └─ Y
│  └─ repeat on (1, 2)
│     └─ Z
├─ +
│  ├─ repeat on (2, 3)
│  │  └─ X
│  ├─ repeat on (2, 3)
│  │  └─ Y
│  └─ repeat on (2, 3)
│     └─ Z
└─ +
   ├─ repeat on (1, 3)
   │  └─ X
   ├─ repeat on (1, 3)
   │  └─ Y
   └─ repeat on (1, 3)
      └─ Z
</pre>


<div class="markdown"><h2>Run a circuit</h2>
</div>


<div class="markdown"><h3>Run circuit on register</h3>
</div>

<pre class='language-julia'><code class='language-julia'>reg = ghz_state(3)</code></pre>
<pre id='var-reg' class='code-output documenter-example-output'>ArrayReg{2, ComplexF64, Array...}
    active qubits: 3/3
    nlevel: 2</pre>

<pre class='language-julia'><code class='language-julia'>@with_terminal print_table(reg)</code></pre>
<pre id="plutouiterminal">
000 ₍₂₎   0.70711 - 0.0im
001 ₍₂₎   0.0 + 0.0im
010 ₍₂₎   0.0 + 0.0im
011 ₍₂₎   0.0 + 0.0im
100 ₍₂₎   0.0 + 0.0im
101 ₍₂₎   0.0 + 0.0im
110 ₍₂₎   0.0 + 0.0im
111 ₍₂₎   0.70711 - 0.0im
</pre>


<pre class='language-julia'><code class='language-julia'>reg2 = apply(reg, qft3)</code></pre>
<pre id='var-reg2' class='code-output documenter-example-output'>ArrayReg{2, ComplexF64, Array...}
    active qubits: 3/3
    nlevel: 2</pre>


<div class="markdown"><p>Note: The inplace version <code>apply&#33;</code> does not allocate and is faster, but Pluto does not track inplace variables very well.</p>
</div>

<pre class='language-julia'><code class='language-julia'>@with_terminal print_table(reg2)</code></pre>
<pre id="plutouiterminal">
000 ₍₂₎   0.5 + 0.0im
001 ₍₂₎   0.42678 - 0.17678im
010 ₍₂₎   0.25 - 0.25im
011 ₍₂₎   0.07322 - 0.17678im
100 ₍₂₎   0.0 + 0.0im
101 ₍₂₎   0.07322 + 0.17678im
110 ₍₂₎   0.25 + 0.25im
111 ₍₂₎   0.42678 + 0.17678im
</pre>



<div class="markdown"><h5>Run a subroutine</h5>
</div>

<pre class='language-julia'><code class='language-julia'>reg16 = rand_state(16)</code></pre>
<pre id='var-reg16' class='code-output documenter-example-output'>ArrayReg{2, ComplexF64, Array...}
    active qubits: 16/16
    nlevel: 2</pre>

<pre class='language-julia'><code class='language-julia'>@time apply(reg16, put(16, 1:8=&gt;EasyBuild.qft_circuit(8)))</code></pre>
<pre id='var-hash172550' class='code-output documenter-example-output'>ArrayReg{2, ComplexF64, Array...}
    active qubits: 16/16
    nlevel: 2</pre>

<pre class='language-julia'><code class='language-julia'>@time apply(reg16, subroutine(16, EasyBuild.qft_circuit(8), 1:8))</code></pre>
<pre id='var-hash437416' class='code-output documenter-example-output'>ArrayReg{2, ComplexF64, Array...}
    active qubits: 16/16
    nlevel: 2</pre>

<pre class='language-julia'><code class='language-julia'>focus!(rand_state(10), (6,2,3)).state</code></pre>
<pre id='var-hash127063' class='code-output documenter-example-output'>8×128 Matrix{ComplexF64}:
 -0.00659563-0.0210159im   -0.0161381+0.00427941im  …  -0.00168342+0.00279972im
  -0.0123904+0.00862976im  -0.0304713-0.0357329im       0.00490122+0.0200389im
  -0.0264094+0.0107274im   -0.0129551-0.0352328im        0.0105061-0.00176262im
  -0.0279175+0.0156832im    0.0104576+0.0146955im        0.0133589-0.00588903im
   0.0192082-0.0236216im   0.00153963+0.0106731im        0.0268599+0.000980262im
 -0.00900309+0.00392223im  -0.0150348+0.040769im    …   -0.0184321-0.00630574im
  -0.0232352-0.0120931im    0.0255792+0.0144375im       -0.0126867-0.0136372im
  -0.0169096+0.0171759im     0.012067-0.0136771im        0.0214318-0.000266297im</pre>


<div class="markdown"><p>A subroutine is similar to <code>put</code>, but runs faster when the sub-circuit size is large. This is because they use different algorithms.</p>
</div>


<div class="markdown"><h3>Run circuit on density matrix</h3>
</div>


<div class="markdown"><p>reduced density can be obtained by calling the <code>density_matrix</code> method.</p>
</div>

<pre class='language-julia'><code class='language-julia'>rho = density_matrix(reg, (2,3))</code></pre>
<pre id='var-rho' class='code-output documenter-example-output'>DensityMatrix{2, ComplexF64, Matrix{ComplexF64}}(ComplexF64[0.5000000000000001 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 - 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 - 0.0im 0.0 - 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 - 0.0im 0.0 - 0.0im 0.0 - 0.0im 0.5000000000000001 + 0.0im])</pre>

<pre class='language-julia'><code class='language-julia'>rho2 = apply(rho, put(2, 2=&gt;X))</code></pre>
<pre id='var-rho2' class='code-output documenter-example-output'>DensityMatrix{2, ComplexF64, Matrix{ComplexF64}}(ComplexF64[0.0 + 0.0im 0.0 + 0.0im 0.0 - 0.0im 0.0 - 0.0im; 0.0 - 0.0im 0.5000000000000001 + 0.0im 0.0 - 0.0im 0.0 - 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.5000000000000001 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 - 0.0im 0.0 + 0.0im])</pre>


<div class="markdown"><p>Density matrix supports unitary channel.</p>
</div>

<pre class='language-julia'><code class='language-julia'>channel = unitary_channel([kron(X,Y), kron(Y,Z), kron(Z,X)], [0.4, 0.3, 0.3])</code></pre>
<pre id='var-channel' class='code-output documenter-example-output'>nqubits: 2
unitary_channel
├─ [0.4] kron
│  ├─ 1=>X
│  └─ 2=>Y
├─ [0.3] kron
│  ├─ 1=>Y
│  └─ 2=>Z
└─ [0.3] kron
   ├─ 1=>Z
   └─ 2=>X
</pre>

<pre class='language-julia'><code class='language-julia'>apply(rho, channel)</code></pre>
<pre id='var-hash430350' class='code-output documenter-example-output'>DensityMatrix{2, ComplexF64, Matrix{ComplexF64}}(ComplexF64[0.20000000000000007 + 0.0im 0.0 + 0.0im 0.0 - 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.30000000000000004 + 0.0im 0.0 + 0.0im 0.0 - 0.0im; 0.0 + 0.0im -0.0 + 0.0im 0.30000000000000004 + 0.0im 0.0 + 0.0im; -0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.20000000000000007 + 0.0im])</pre>


<div class="markdown"><p>Here, <code>kron</code> is a composite block that defined as</p>
<p class="tex">$$\texttt&#123;kron&#125;&#40;G_1, G_2, \ldots, G_n&#41; &#61; G_n \otimes G_&#123;n-1&#125; \otimes \ldots \otimes G_1$$</p>
<p>A <code>unitary_channel</code> does not have a matrix representation. When applying a unitary channel <code>unitary_channel&#40;U, p&#41;</code> to a density matrix, it effectively does</p>
<p class="tex">$$\phi&#40;\rho&#41; &#61; \sum_i p_i U_i ρ U_i^\dagger$$</p>
</div>


<div class="markdown"><h2>Manipulate the block tree</h2>
</div>


<div class="markdown"><h3>Recurse over the block tree</h3>
</div>


<div class="markdown"><p>Yao&#39;s block representation &#40;QBIR&#41; defines a tree, its siblings can be accessed with the <code>subblocks</code> function.</p>
</div>

<pre class='language-julia'><code class='language-julia'>subblocks(qft3)</code></pre>
<pre id='var-hash572067' class='code-output documenter-example-output'>3-element Vector{AbstractBlock{2}}:
 nqubits: 3
chain
├─ put on (1)
│  └─ H
├─ control(2)
│  └─ (1,) shift(1.5707963267948966)
└─ control(3)
   └─ (1,) shift(0.7853981633974483)

 nqubits: 3
chain
├─ put on (2)
│  └─ H
└─ control(3)
   └─ (2,) shift(1.5707963267948966)

 nqubits: 3
chain
└─ put on (3)
   └─ H
</pre>

<pre class='language-julia'><code class='language-julia'>subblocks(subblocks(qft3)[1])</code></pre>
<pre id='var-hash146294' class='code-output documenter-example-output'>3-element Vector{AbstractBlock{2}}:
 nqubits: 3
put on (1)
└─ H
 nqubits: 3
control(2)
└─ (1,) shift(1.5707963267948966)
 nqubits: 3
control(3)
└─ (1,) shift(0.7853981633974483)</pre>

<pre class='language-julia'><code class='language-julia'>subblocks(subblocks(subblocks(qft3)[1])[1])</code></pre>
<pre id='var-hash411747' class='code-output documenter-example-output'>(H,)</pre>

<pre class='language-julia'><code class='language-julia'>subblocks(subblocks(subblocks(subblocks(qft3)[1])[1])[1])</code></pre>
<pre id='var-hash127716' class='code-output documenter-example-output'>()</pre>

<pre class='language-julia'><code class='language-julia'>begin
    decompose(x::HGate)= Rz(0.5π)*Rx(0.5π)*Rz(0.5π)
    decompose(x::AbstractBlock)= chsubblocks(x, decompose.(subblocks(x)))
    vizcircuit(decompose(qft3); scale=0.65)
end</code></pre>
<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:xlink="http://www.w3.org/1999/xlink"
     version="1.2"
     width="179.56mm" height="36.56mm" viewBox="0 0 179.56 36.56"
     stroke="none"
     fill="#000000"
     stroke-width="0.3"
     font-size="3.88"

     id="img-39279afb">
<defs>
  <marker id="arrow" markerWidth="15" markerHeight="7" refX="5" refY="3.5" orient="auto" markerUnits="strokeWidth">
    <path d="M0,0 L15,3.5 L0,7 z" stroke="context-stroke" fill="context-stroke"/>
  </marker>
</defs>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-39279afb-1">
  <g transform="translate(166.97,28.44)">
    <path d="M-5.4,-3.25 L5.4,-3.25 5.4,3.25 -5.4,3.25  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-39279afb-2">
  <g transform="translate(153.16,28.44)">
    <path d="M-5.4,-3.25 L5.4,-3.25 5.4,3.25 -5.4,3.25  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-39279afb-3">
  <g transform="translate(139.34,28.44)">
    <path d="M-5.4,-3.25 L5.4,-3.25 5.4,3.25 -5.4,3.25  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0" fill="#000000" id="img-39279afb-4">
  <g transform="translate(125.53,28.44)">
    <circle cx="0" cy="0" r="0.54" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-39279afb-5">
  <g transform="translate(125.53,16.25)">
    <path d="M-4.68,-3.25 L4.68,-3.25 4.68,3.25 -4.68,3.25  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-39279afb-6">
  <g transform="translate(111.72,16.25)">
    <path d="M-5.4,-3.25 L5.4,-3.25 5.4,3.25 -5.4,3.25  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-39279afb-7">
  <g transform="translate(97.91,16.25)">
    <path d="M-5.4,-3.25 L5.4,-3.25 5.4,3.25 -5.4,3.25  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-39279afb-8">
  <g transform="translate(84.09,16.25)">
    <path d="M-5.4,-3.25 L5.4,-3.25 5.4,3.25 -5.4,3.25  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0" fill="#000000" id="img-39279afb-9">
  <g transform="translate(70.28,28.44)">
    <circle cx="0" cy="0" r="0.54" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-39279afb-10">
  <g transform="translate(70.28,4.06)">
    <path d="M-4.68,-3.25 L4.68,-3.25 4.68,3.25 -4.68,3.25  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0" fill="#000000" id="img-39279afb-11">
  <g transform="translate(56.47,16.25)">
    <circle cx="0" cy="0" r="0.54" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-39279afb-12">
  <g transform="translate(56.47,4.06)">
    <path d="M-4.68,-3.25 L4.68,-3.25 4.68,3.25 -4.68,3.25  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-39279afb-13">
  <g transform="translate(42.66,4.06)">
    <path d="M-5.4,-3.25 L5.4,-3.25 5.4,3.25 -5.4,3.25  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-39279afb-14">
  <g transform="translate(28.84,4.06)">
    <path d="M-5.4,-3.25 L5.4,-3.25 5.4,3.25 -5.4,3.25  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" fill="#000000" fill-opacity="0.000" id="img-39279afb-15">
  <g transform="translate(15.03,4.06)">
    <path d="M-5.4,-3.25 L5.4,-3.25 5.4,3.25 -5.4,3.25  z" class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-16">
  <g transform="translate(173.88,28.44)">
    <path fill="none" d="M0,0 L0,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-17">
  <g transform="translate(153.16,16.25)">
    <path fill="none" d="M20.72,0 L-20.72,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-18">
  <g transform="translate(125.53,4.06)">
    <path fill="none" d="M48.34,0 L-48.34,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-19">
  <g transform="translate(173.12,28.44)">
    <path fill="none" d="M-0.75,0 L0.75,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-20">
  <g transform="translate(160.82,28.44)">
    <path fill="none" d="M-0.75,0 L0.75,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-21">
  <g transform="translate(159.31,28.44)">
    <path fill="none" d="M-0.75,0 L0.75,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-22">
  <g transform="translate(147,28.44)">
    <path fill="none" d="M-0.75,0 L0.75,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-23">
  <g transform="translate(145.5,28.44)">
    <path fill="none" d="M-0.75,0 L0.75,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-24">
  <g transform="translate(133.19,28.44)">
    <path fill="none" d="M-0.75,0 L0.75,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-25">
  <g transform="translate(129.26,28.44)">
    <path fill="none" d="M-3.18,0 L3.18,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-26">
  <g transform="translate(101.09,28.44)">
    <path fill="none" d="M-23.9,0 L23.9,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-27">
  <g transform="translate(131.33,16.25)">
    <path fill="none" d="M-1.11,0 L1.11,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-28">
  <g transform="translate(119.74,16.25)">
    <path fill="none" d="M-1.11,0 L1.11,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-29">
  <g transform="translate(125.53,23.7)">
    <path fill="none" d="M0,4.2 L0,-4.2 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-30">
  <g transform="translate(117.87,16.25)">
    <path fill="none" d="M-0.75,0 L0.75,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-31">
  <g transform="translate(105.57,16.25)">
    <path fill="none" d="M-0.75,0 L0.75,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-32">
  <g transform="translate(104.06,16.25)">
    <path fill="none" d="M-0.75,0 L0.75,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-33">
  <g transform="translate(91.75,16.25)">
    <path fill="none" d="M-0.75,0 L0.75,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-34">
  <g transform="translate(90.25,16.25)">
    <path fill="none" d="M-0.75,0 L0.75,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-35">
  <g transform="translate(77.94,16.25)">
    <path fill="none" d="M-0.75,0 L0.75,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-36">
  <g transform="translate(70.28,16.25)">
    <path fill="none" d="M-6.91,0 L6.91,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-37">
  <g transform="translate(74.01,28.44)">
    <path fill="none" d="M-3.18,0 L3.18,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-38">
  <g transform="translate(38.93,28.44)">
    <path fill="none" d="M-30.81,0 L30.81,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-39">
  <g transform="translate(76.08,4.06)">
    <path fill="none" d="M-1.11,0 L1.11,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-40">
  <g transform="translate(64.49,4.06)">
    <path fill="none" d="M-1.11,0 L1.11,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-41">
  <g transform="translate(70.28,17.6)">
    <path fill="none" d="M0,10.29 L0,-10.29 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-42">
  <g transform="translate(60.19,16.25)">
    <path fill="none" d="M-3.18,0 L3.18,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-43">
  <g transform="translate(32.03,16.25)">
    <path fill="none" d="M-23.9,0 L23.9,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-44">
  <g transform="translate(62.26,4.06)">
    <path fill="none" d="M-1.11,0 L1.11,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-45">
  <g transform="translate(50.67,4.06)">
    <path fill="none" d="M-1.11,0 L1.11,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-46">
  <g transform="translate(56.47,11.51)">
    <path fill="none" d="M0,4.2 L0,-4.2 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-47">
  <g transform="translate(48.81,4.06)">
    <path fill="none" d="M-0.75,0 L0.75,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-48">
  <g transform="translate(36.5,4.06)">
    <path fill="none" d="M-0.75,0 L0.75,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-49">
  <g transform="translate(35,4.06)">
    <path fill="none" d="M-0.75,0 L0.75,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-50">
  <g transform="translate(22.69,4.06)">
    <path fill="none" d="M-0.75,0 L0.75,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-51">
  <g transform="translate(21.18,4.06)">
    <path fill="none" d="M-0.75,0 L0.75,0 " class="primitive"/>
  </g>
</g>
<g stroke-width="0.35" stroke="#000000" id="img-39279afb-52">
  <g transform="translate(8.88,4.06)">
    <path fill="none" d="M-0.75,0 L0.75,0 " class="primitive"/>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="2.29" id="img-39279afb-53">
  <g transform="translate(166.97,28.44)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">Rz(π/2)</text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="2.29" id="img-39279afb-54">
  <g transform="translate(153.16,28.44)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">Rx(π/2)</text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="2.29" id="img-39279afb-55">
  <g transform="translate(139.34,28.44)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">Rz(π/2)</text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="3.67" id="img-39279afb-56">
  <g transform="translate(125.53,28.44)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em"></text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="2.29" id="img-39279afb-57">
  <g transform="translate(125.53,16.25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">ϕ(π/2)</text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="2.29" id="img-39279afb-58">
  <g transform="translate(111.72,16.25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">Rz(π/2)</text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="2.29" id="img-39279afb-59">
  <g transform="translate(97.91,16.25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">Rx(π/2)</text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="2.29" id="img-39279afb-60">
  <g transform="translate(84.09,16.25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">Rz(π/2)</text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="3.67" id="img-39279afb-61">
  <g transform="translate(70.28,28.44)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em"></text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="2.29" id="img-39279afb-62">
  <g transform="translate(70.28,4.06)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">ϕ(π/4)</text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="3.67" id="img-39279afb-63">
  <g transform="translate(56.47,16.25)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em"></text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="2.29" id="img-39279afb-64">
  <g transform="translate(56.47,4.06)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">ϕ(π/2)</text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="2.29" id="img-39279afb-65">
  <g transform="translate(42.66,4.06)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">Rz(π/2)</text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="2.29" id="img-39279afb-66">
  <g transform="translate(28.84,4.06)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">Rx(π/2)</text>
    </g>
  </g>
</g>
<g font-family="monospace" fill="#000000" font-size="2.29" id="img-39279afb-67">
  <g transform="translate(15.03,4.06)">
    <g class="primitive">
      <text text-anchor="middle" dy="0.35em">Rz(π/2)</text>
    </g>
  </g>
</g>
<script> <![CDATA[
(function(N){var k=/[\.\/]/,L=/\s*,\s*/,C=function(a,d){return a-d},a,v,y={n:{}},M=function(){for(var a=0,d=this.length;a<d;a++)if("undefined"!=typeof this[a])return this[a]},A=function(){for(var a=this.length;--a;)if("undefined"!=typeof this[a])return this[a]},w=function(k,d){k=String(k);var f=v,n=Array.prototype.slice.call(arguments,2),u=w.listeners(k),p=0,b,q=[],e={},l=[],r=a;l.firstDefined=M;l.lastDefined=A;a=k;for(var s=v=0,x=u.length;s<x;s++)"zIndex"in u[s]&&(q.push(u[s].zIndex),0>u[s].zIndex&&
(e[u[s].zIndex]=u[s]));for(q.sort(C);0>q[p];)if(b=e[q[p++] ],l.push(b.apply(d,n)),v)return v=f,l;for(s=0;s<x;s++)if(b=u[s],"zIndex"in b)if(b.zIndex==q[p]){l.push(b.apply(d,n));if(v)break;do if(p++,(b=e[q[p] ])&&l.push(b.apply(d,n)),v)break;while(b)}else e[b.zIndex]=b;else if(l.push(b.apply(d,n)),v)break;v=f;a=r;return l};w._events=y;w.listeners=function(a){a=a.split(k);var d=y,f,n,u,p,b,q,e,l=[d],r=[];u=0;for(p=a.length;u<p;u++){e=[];b=0;for(q=l.length;b<q;b++)for(d=l[b].n,f=[d[a[u] ],d["*"] ],n=2;n--;)if(d=
f[n])e.push(d),r=r.concat(d.f||[]);l=e}return r};w.on=function(a,d){a=String(a);if("function"!=typeof d)return function(){};for(var f=a.split(L),n=0,u=f.length;n<u;n++)(function(a){a=a.split(k);for(var b=y,f,e=0,l=a.length;e<l;e++)b=b.n,b=b.hasOwnProperty(a[e])&&b[a[e] ]||(b[a[e] ]={n:{}});b.f=b.f||[];e=0;for(l=b.f.length;e<l;e++)if(b.f[e]==d){f=!0;break}!f&&b.f.push(d)})(f[n]);return function(a){+a==+a&&(d.zIndex=+a)}};w.f=function(a){var d=[].slice.call(arguments,1);return function(){w.apply(null,
[a,null].concat(d).concat([].slice.call(arguments,0)))}};w.stop=function(){v=1};w.nt=function(k){return k?(new RegExp("(?:\\.|\\/|^)"+k+"(?:\\.|\\/|$)")).test(a):a};w.nts=function(){return a.split(k)};w.off=w.unbind=function(a,d){if(a){var f=a.split(L);if(1<f.length)for(var n=0,u=f.length;n<u;n++)w.off(f[n],d);else{for(var f=a.split(k),p,b,q,e,l=[y],n=0,u=f.length;n<u;n++)for(e=0;e<l.length;e+=q.length-2){q=[e,1];p=l[e].n;if("*"!=f[n])p[f[n] ]&&q.push(p[f[n] ]);else for(b in p)p.hasOwnProperty(b)&&
q.push(p[b]);l.splice.apply(l,q)}n=0;for(u=l.length;n<u;n++)for(p=l[n];p.n;){if(d){if(p.f){e=0;for(f=p.f.length;e<f;e++)if(p.f[e]==d){p.f.splice(e,1);break}!p.f.length&&delete p.f}for(b in p.n)if(p.n.hasOwnProperty(b)&&p.n[b].f){q=p.n[b].f;e=0;for(f=q.length;e<f;e++)if(q[e]==d){q.splice(e,1);break}!q.length&&delete p.n[b].f}}else for(b in delete p.f,p.n)p.n.hasOwnProperty(b)&&p.n[b].f&&delete p.n[b].f;p=p.n}}}else w._events=y={n:{}}};w.once=function(a,d){var f=function(){w.unbind(a,f);return d.apply(this,
arguments)};return w.on(a,f)};w.version="0.4.2";w.toString=function(){return"You are running Eve 0.4.2"};"undefined"!=typeof module&&module.exports?module.exports=w:"function"===typeof define&&define.amd?define("eve",[],function(){return w}):N.eve=w})(this);
(function(N,k){"function"===typeof define&&define.amd?define("Snap.svg",["eve"],function(L){return k(N,L)}):k(N,N.eve)})(this,function(N,k){var L=function(a){var k={},y=N.requestAnimationFrame||N.webkitRequestAnimationFrame||N.mozRequestAnimationFrame||N.oRequestAnimationFrame||N.msRequestAnimationFrame||function(a){setTimeout(a,16)},M=Array.isArray||function(a){return a instanceof Array||"[object Array]"==Object.prototype.toString.call(a)},A=0,w="M"+(+new Date).toString(36),z=function(a){if(null==
a)return this.s;var b=this.s-a;this.b+=this.dur*b;this.B+=this.dur*b;this.s=a},d=function(a){if(null==a)return this.spd;this.spd=a},f=function(a){if(null==a)return this.dur;this.s=this.s*a/this.dur;this.dur=a},n=function(){delete k[this.id];this.update();a("mina.stop."+this.id,this)},u=function(){this.pdif||(delete k[this.id],this.update(),this.pdif=this.get()-this.b)},p=function(){this.pdif&&(this.b=this.get()-this.pdif,delete this.pdif,k[this.id]=this)},b=function(){var a;if(M(this.start)){a=[];
for(var b=0,e=this.start.length;b<e;b++)a[b]=+this.start[b]+(this.end[b]-this.start[b])*this.easing(this.s)}else a=+this.start+(this.end-this.start)*this.easing(this.s);this.set(a)},q=function(){var l=0,b;for(b in k)if(k.hasOwnProperty(b)){var e=k[b],f=e.get();l++;e.s=(f-e.b)/(e.dur/e.spd);1<=e.s&&(delete k[b],e.s=1,l--,function(b){setTimeout(function(){a("mina.finish."+b.id,b)})}(e));e.update()}l&&y(q)},e=function(a,r,s,x,G,h,J){a={id:w+(A++).toString(36),start:a,end:r,b:s,s:0,dur:x-s,spd:1,get:G,
set:h,easing:J||e.linear,status:z,speed:d,duration:f,stop:n,pause:u,resume:p,update:b};k[a.id]=a;r=0;for(var K in k)if(k.hasOwnProperty(K)&&(r++,2==r))break;1==r&&y(q);return a};e.time=Date.now||function(){return+new Date};e.getById=function(a){return k[a]||null};e.linear=function(a){return a};e.easeout=function(a){return Math.pow(a,1.7)};e.easein=function(a){return Math.pow(a,0.48)};e.easeinout=function(a){if(1==a)return 1;if(0==a)return 0;var b=0.48-a/1.04,e=Math.sqrt(0.1734+b*b);a=e-b;a=Math.pow(Math.abs(a),
1/3)*(0>a?-1:1);b=-e-b;b=Math.pow(Math.abs(b),1/3)*(0>b?-1:1);a=a+b+0.5;return 3*(1-a)*a*a+a*a*a};e.backin=function(a){return 1==a?1:a*a*(2.70158*a-1.70158)};e.backout=function(a){if(0==a)return 0;a-=1;return a*a*(2.70158*a+1.70158)+1};e.elastic=function(a){return a==!!a?a:Math.pow(2,-10*a)*Math.sin(2*(a-0.075)*Math.PI/0.3)+1};e.bounce=function(a){a<1/2.75?a*=7.5625*a:a<2/2.75?(a-=1.5/2.75,a=7.5625*a*a+0.75):a<2.5/2.75?(a-=2.25/2.75,a=7.5625*a*a+0.9375):(a-=2.625/2.75,a=7.5625*a*a+0.984375);return a};
return N.mina=e}("undefined"==typeof k?function(){}:k),C=function(){function a(c,t){if(c){if(c.tagName)return x(c);if(y(c,"array")&&a.set)return a.set.apply(a,c);if(c instanceof e)return c;if(null==t)return c=G.doc.querySelector(c),x(c)}return new s(null==c?"100%":c,null==t?"100%":t)}function v(c,a){if(a){"#text"==c&&(c=G.doc.createTextNode(a.text||""));"string"==typeof c&&(c=v(c));if("string"==typeof a)return"xlink:"==a.substring(0,6)?c.getAttributeNS(m,a.substring(6)):"xml:"==a.substring(0,4)?c.getAttributeNS(la,
a.substring(4)):c.getAttribute(a);for(var da in a)if(a[h](da)){var b=J(a[da]);b?"xlink:"==da.substring(0,6)?c.setAttributeNS(m,da.substring(6),b):"xml:"==da.substring(0,4)?c.setAttributeNS(la,da.substring(4),b):c.setAttribute(da,b):c.removeAttribute(da)}}else c=G.doc.createElementNS(la,c);return c}function y(c,a){a=J.prototype.toLowerCase.call(a);return"finite"==a?isFinite(c):"array"==a&&(c instanceof Array||Array.isArray&&Array.isArray(c))?!0:"null"==a&&null===c||a==typeof c&&null!==c||"object"==
a&&c===Object(c)||$.call(c).slice(8,-1).toLowerCase()==a}function M(c){if("function"==typeof c||Object(c)!==c)return c;var a=new c.constructor,b;for(b in c)c[h](b)&&(a[b]=M(c[b]));return a}function A(c,a,b){function m(){var e=Array.prototype.slice.call(arguments,0),f=e.join("\u2400"),d=m.cache=m.cache||{},l=m.count=m.count||[];if(d[h](f)){a:for(var e=l,l=f,B=0,H=e.length;B<H;B++)if(e[B]===l){e.push(e.splice(B,1)[0]);break a}return b?b(d[f]):d[f]}1E3<=l.length&&delete d[l.shift()];l.push(f);d[f]=c.apply(a,
e);return b?b(d[f]):d[f]}return m}function w(c,a,b,m,e,f){return null==e?(c-=b,a-=m,c||a?(180*I.atan2(-a,-c)/C+540)%360:0):w(c,a,e,f)-w(b,m,e,f)}function z(c){return c%360*C/180}function d(c){var a=[];c=c.replace(/(?:^|\s)(\w+)\(([^)]+)\)/g,function(c,b,m){m=m.split(/\s*,\s*|\s+/);"rotate"==b&&1==m.length&&m.push(0,0);"scale"==b&&(2<m.length?m=m.slice(0,2):2==m.length&&m.push(0,0),1==m.length&&m.push(m[0],0,0));"skewX"==b?a.push(["m",1,0,I.tan(z(m[0])),1,0,0]):"skewY"==b?a.push(["m",1,I.tan(z(m[0])),
0,1,0,0]):a.push([b.charAt(0)].concat(m));return c});return a}function f(c,t){var b=O(c),m=new a.Matrix;if(b)for(var e=0,f=b.length;e<f;e++){var h=b[e],d=h.length,B=J(h[0]).toLowerCase(),H=h[0]!=B,l=H?m.invert():0,E;"t"==B&&2==d?m.translate(h[1],0):"t"==B&&3==d?H?(d=l.x(0,0),B=l.y(0,0),H=l.x(h[1],h[2]),l=l.y(h[1],h[2]),m.translate(H-d,l-B)):m.translate(h[1],h[2]):"r"==B?2==d?(E=E||t,m.rotate(h[1],E.x+E.width/2,E.y+E.height/2)):4==d&&(H?(H=l.x(h[2],h[3]),l=l.y(h[2],h[3]),m.rotate(h[1],H,l)):m.rotate(h[1],
h[2],h[3])):"s"==B?2==d||3==d?(E=E||t,m.scale(h[1],h[d-1],E.x+E.width/2,E.y+E.height/2)):4==d?H?(H=l.x(h[2],h[3]),l=l.y(h[2],h[3]),m.scale(h[1],h[1],H,l)):m.scale(h[1],h[1],h[2],h[3]):5==d&&(H?(H=l.x(h[3],h[4]),l=l.y(h[3],h[4]),m.scale(h[1],h[2],H,l)):m.scale(h[1],h[2],h[3],h[4])):"m"==B&&7==d&&m.add(h[1],h[2],h[3],h[4],h[5],h[6])}return m}function n(c,t){if(null==t){var m=!0;t="linearGradient"==c.type||"radialGradient"==c.type?c.node.getAttribute("gradientTransform"):"pattern"==c.type?c.node.getAttribute("patternTransform"):
c.node.getAttribute("transform");if(!t)return new a.Matrix;t=d(t)}else t=a._.rgTransform.test(t)?J(t).replace(/\.{3}|\u2026/g,c._.transform||aa):d(t),y(t,"array")&&(t=a.path?a.path.toString.call(t):J(t)),c._.transform=t;var b=f(t,c.getBBox(1));if(m)return b;c.matrix=b}function u(c){c=c.node.ownerSVGElement&&x(c.node.ownerSVGElement)||c.node.parentNode&&x(c.node.parentNode)||a.select("svg")||a(0,0);var t=c.select("defs"),t=null==t?!1:t.node;t||(t=r("defs",c.node).node);return t}function p(c){return c.node.ownerSVGElement&&
x(c.node.ownerSVGElement)||a.select("svg")}function b(c,a,m){function b(c){if(null==c)return aa;if(c==+c)return c;v(B,{width:c});try{return B.getBBox().width}catch(a){return 0}}function h(c){if(null==c)return aa;if(c==+c)return c;v(B,{height:c});try{return B.getBBox().height}catch(a){return 0}}function e(b,B){null==a?d[b]=B(c.attr(b)||0):b==a&&(d=B(null==m?c.attr(b)||0:m))}var f=p(c).node,d={},B=f.querySelector(".svg---mgr");B||(B=v("rect"),v(B,{x:-9E9,y:-9E9,width:10,height:10,"class":"svg---mgr",
fill:"none"}),f.appendChild(B));switch(c.type){case "rect":e("rx",b),e("ry",h);case "image":e("width",b),e("height",h);case "text":e("x",b);e("y",h);break;case "circle":e("cx",b);e("cy",h);e("r",b);break;case "ellipse":e("cx",b);e("cy",h);e("rx",b);e("ry",h);break;case "line":e("x1",b);e("x2",b);e("y1",h);e("y2",h);break;case "marker":e("refX",b);e("markerWidth",b);e("refY",h);e("markerHeight",h);break;case "radialGradient":e("fx",b);e("fy",h);break;case "tspan":e("dx",b);e("dy",h);break;default:e(a,
b)}f.removeChild(B);return d}function q(c){y(c,"array")||(c=Array.prototype.slice.call(arguments,0));for(var a=0,b=0,m=this.node;this[a];)delete this[a++];for(a=0;a<c.length;a++)"set"==c[a].type?c[a].forEach(function(c){m.appendChild(c.node)}):m.appendChild(c[a].node);for(var h=m.childNodes,a=0;a<h.length;a++)this[b++]=x(h[a]);return this}function e(c){if(c.snap in E)return E[c.snap];var a=this.id=V(),b;try{b=c.ownerSVGElement}catch(m){}this.node=c;b&&(this.paper=new s(b));this.type=c.tagName;this.anims=
{};this._={transform:[]};c.snap=a;E[a]=this;"g"==this.type&&(this.add=q);if(this.type in{g:1,mask:1,pattern:1})for(var e in s.prototype)s.prototype[h](e)&&(this[e]=s.prototype[e])}function l(c){this.node=c}function r(c,a){var b=v(c);a.appendChild(b);return x(b)}function s(c,a){var b,m,f,d=s.prototype;if(c&&"svg"==c.tagName){if(c.snap in E)return E[c.snap];var l=c.ownerDocument;b=new e(c);m=c.getElementsByTagName("desc")[0];f=c.getElementsByTagName("defs")[0];m||(m=v("desc"),m.appendChild(l.createTextNode("Created with Snap")),
b.node.appendChild(m));f||(f=v("defs"),b.node.appendChild(f));b.defs=f;for(var ca in d)d[h](ca)&&(b[ca]=d[ca]);b.paper=b.root=b}else b=r("svg",G.doc.body),v(b.node,{height:a,version:1.1,width:c,xmlns:la});return b}function x(c){return!c||c instanceof e||c instanceof l?c:c.tagName&&"svg"==c.tagName.toLowerCase()?new s(c):c.tagName&&"object"==c.tagName.toLowerCase()&&"image/svg+xml"==c.type?new s(c.contentDocument.getElementsByTagName("svg")[0]):new e(c)}a.version="0.3.0";a.toString=function(){return"Snap v"+
this.version};a._={};var G={win:N,doc:N.document};a._.glob=G;var h="hasOwnProperty",J=String,K=parseFloat,U=parseInt,I=Math,P=I.max,Q=I.min,Y=I.abs,C=I.PI,aa="",$=Object.prototype.toString,F=/^\s*((#[a-f\d]{6})|(#[a-f\d]{3})|rgba?\(\s*([\d\.]+%?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+%?(?:\s*,\s*[\d\.]+%?)?)\s*\)|hsba?\(\s*([\d\.]+(?:deg|\xb0|%)?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+(?:%?\s*,\s*[\d\.]+)?%?)\s*\)|hsla?\(\s*([\d\.]+(?:deg|\xb0|%)?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+(?:%?\s*,\s*[\d\.]+)?%?)\s*\))\s*$/i;a._.separator=
RegExp("[,\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]+");var S=RegExp("[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*"),X={hs:1,rg:1},W=RegExp("([a-z])[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029,]*((-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*)+)",
"ig"),ma=RegExp("([rstm])[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029,]*((-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*)+)","ig"),Z=RegExp("(-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?)[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*",
"ig"),na=0,ba="S"+(+new Date).toString(36),V=function(){return ba+(na++).toString(36)},m="http://www.w3.org/1999/xlink",la="http://www.w3.org/2000/svg",E={},ca=a.url=function(c){return"url('#"+c+"')"};a._.$=v;a._.id=V;a.format=function(){var c=/\{([^\}]+)\}/g,a=/(?:(?:^|\.)(.+?)(?=\[|\.|$|\()|\[('|")(.+?)\2\])(\(\))?/g,b=function(c,b,m){var h=m;b.replace(a,function(c,a,b,m,t){a=a||m;h&&(a in h&&(h=h[a]),"function"==typeof h&&t&&(h=h()))});return h=(null==h||h==m?c:h)+""};return function(a,m){return J(a).replace(c,
function(c,a){return b(c,a,m)})}}();a._.clone=M;a._.cacher=A;a.rad=z;a.deg=function(c){return 180*c/C%360};a.angle=w;a.is=y;a.snapTo=function(c,a,b){b=y(b,"finite")?b:10;if(y(c,"array"))for(var m=c.length;m--;){if(Y(c[m]-a)<=b)return c[m]}else{c=+c;m=a%c;if(m<b)return a-m;if(m>c-b)return a-m+c}return a};a.getRGB=A(function(c){if(!c||(c=J(c)).indexOf("-")+1)return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka};if("none"==c)return{r:-1,g:-1,b:-1,hex:"none",toString:ka};!X[h](c.toLowerCase().substring(0,
2))&&"#"!=c.charAt()&&(c=T(c));if(!c)return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka};var b,m,e,f,d;if(c=c.match(F)){c[2]&&(e=U(c[2].substring(5),16),m=U(c[2].substring(3,5),16),b=U(c[2].substring(1,3),16));c[3]&&(e=U((d=c[3].charAt(3))+d,16),m=U((d=c[3].charAt(2))+d,16),b=U((d=c[3].charAt(1))+d,16));c[4]&&(d=c[4].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b*=2.55),m=K(d[1]),"%"==d[1].slice(-1)&&(m*=2.55),e=K(d[2]),"%"==d[2].slice(-1)&&(e*=2.55),"rgba"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),
d[3]&&"%"==d[3].slice(-1)&&(f/=100));if(c[5])return d=c[5].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b/=100),m=K(d[1]),"%"==d[1].slice(-1)&&(m/=100),e=K(d[2]),"%"==d[2].slice(-1)&&(e/=100),"deg"!=d[0].slice(-3)&&"\u00b0"!=d[0].slice(-1)||(b/=360),"hsba"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),d[3]&&"%"==d[3].slice(-1)&&(f/=100),a.hsb2rgb(b,m,e,f);if(c[6])return d=c[6].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b/=100),m=K(d[1]),"%"==d[1].slice(-1)&&(m/=100),e=K(d[2]),"%"==d[2].slice(-1)&&(e/=100),
"deg"!=d[0].slice(-3)&&"\u00b0"!=d[0].slice(-1)||(b/=360),"hsla"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),d[3]&&"%"==d[3].slice(-1)&&(f/=100),a.hsl2rgb(b,m,e,f);b=Q(I.round(b),255);m=Q(I.round(m),255);e=Q(I.round(e),255);f=Q(P(f,0),1);c={r:b,g:m,b:e,toString:ka};c.hex="#"+(16777216|e|m<<8|b<<16).toString(16).slice(1);c.opacity=y(f,"finite")?f:1;return c}return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka}},a);a.hsb=A(function(c,b,m){return a.hsb2rgb(c,b,m).hex});a.hsl=A(function(c,b,m){return a.hsl2rgb(c,
b,m).hex});a.rgb=A(function(c,a,b,m){if(y(m,"finite")){var e=I.round;return"rgba("+[e(c),e(a),e(b),+m.toFixed(2)]+")"}return"#"+(16777216|b|a<<8|c<<16).toString(16).slice(1)});var T=function(c){var a=G.doc.getElementsByTagName("head")[0]||G.doc.getElementsByTagName("svg")[0];T=A(function(c){if("red"==c.toLowerCase())return"rgb(255, 0, 0)";a.style.color="rgb(255, 0, 0)";a.style.color=c;c=G.doc.defaultView.getComputedStyle(a,aa).getPropertyValue("color");return"rgb(255, 0, 0)"==c?null:c});return T(c)},
qa=function(){return"hsb("+[this.h,this.s,this.b]+")"},ra=function(){return"hsl("+[this.h,this.s,this.l]+")"},ka=function(){return 1==this.opacity||null==this.opacity?this.hex:"rgba("+[this.r,this.g,this.b,this.opacity]+")"},D=function(c,b,m){null==b&&y(c,"object")&&"r"in c&&"g"in c&&"b"in c&&(m=c.b,b=c.g,c=c.r);null==b&&y(c,string)&&(m=a.getRGB(c),c=m.r,b=m.g,m=m.b);if(1<c||1<b||1<m)c/=255,b/=255,m/=255;return[c,b,m]},oa=function(c,b,m,e){c=I.round(255*c);b=I.round(255*b);m=I.round(255*m);c={r:c,
g:b,b:m,opacity:y(e,"finite")?e:1,hex:a.rgb(c,b,m),toString:ka};y(e,"finite")&&(c.opacity=e);return c};a.color=function(c){var b;y(c,"object")&&"h"in c&&"s"in c&&"b"in c?(b=a.hsb2rgb(c),c.r=b.r,c.g=b.g,c.b=b.b,c.opacity=1,c.hex=b.hex):y(c,"object")&&"h"in c&&"s"in c&&"l"in c?(b=a.hsl2rgb(c),c.r=b.r,c.g=b.g,c.b=b.b,c.opacity=1,c.hex=b.hex):(y(c,"string")&&(c=a.getRGB(c)),y(c,"object")&&"r"in c&&"g"in c&&"b"in c&&!("error"in c)?(b=a.rgb2hsl(c),c.h=b.h,c.s=b.s,c.l=b.l,b=a.rgb2hsb(c),c.v=b.b):(c={hex:"none"},
c.r=c.g=c.b=c.h=c.s=c.v=c.l=-1,c.error=1));c.toString=ka;return c};a.hsb2rgb=function(c,a,b,m){y(c,"object")&&"h"in c&&"s"in c&&"b"in c&&(b=c.b,a=c.s,c=c.h,m=c.o);var e,h,d;c=360*c%360/60;d=b*a;a=d*(1-Y(c%2-1));b=e=h=b-d;c=~~c;b+=[d,a,0,0,a,d][c];e+=[a,d,d,a,0,0][c];h+=[0,0,a,d,d,a][c];return oa(b,e,h,m)};a.hsl2rgb=function(c,a,b,m){y(c,"object")&&"h"in c&&"s"in c&&"l"in c&&(b=c.l,a=c.s,c=c.h);if(1<c||1<a||1<b)c/=360,a/=100,b/=100;var e,h,d;c=360*c%360/60;d=2*a*(0.5>b?b:1-b);a=d*(1-Y(c%2-1));b=e=
h=b-d/2;c=~~c;b+=[d,a,0,0,a,d][c];e+=[a,d,d,a,0,0][c];h+=[0,0,a,d,d,a][c];return oa(b,e,h,m)};a.rgb2hsb=function(c,a,b){b=D(c,a,b);c=b[0];a=b[1];b=b[2];var m,e;m=P(c,a,b);e=m-Q(c,a,b);c=((0==e?0:m==c?(a-b)/e:m==a?(b-c)/e+2:(c-a)/e+4)+360)%6*60/360;return{h:c,s:0==e?0:e/m,b:m,toString:qa}};a.rgb2hsl=function(c,a,b){b=D(c,a,b);c=b[0];a=b[1];b=b[2];var m,e,h;m=P(c,a,b);e=Q(c,a,b);h=m-e;c=((0==h?0:m==c?(a-b)/h:m==a?(b-c)/h+2:(c-a)/h+4)+360)%6*60/360;m=(m+e)/2;return{h:c,s:0==h?0:0.5>m?h/(2*m):h/(2-2*
m),l:m,toString:ra}};a.parsePathString=function(c){if(!c)return null;var b=a.path(c);if(b.arr)return a.path.clone(b.arr);var m={a:7,c:6,o:2,h:1,l:2,m:2,r:4,q:4,s:4,t:2,v:1,u:3,z:0},e=[];y(c,"array")&&y(c[0],"array")&&(e=a.path.clone(c));e.length||J(c).replace(W,function(c,a,b){var h=[];c=a.toLowerCase();b.replace(Z,function(c,a){a&&h.push(+a)});"m"==c&&2<h.length&&(e.push([a].concat(h.splice(0,2))),c="l",a="m"==a?"l":"L");"o"==c&&1==h.length&&e.push([a,h[0] ]);if("r"==c)e.push([a].concat(h));else for(;h.length>=
m[c]&&(e.push([a].concat(h.splice(0,m[c]))),m[c]););});e.toString=a.path.toString;b.arr=a.path.clone(e);return e};var O=a.parseTransformString=function(c){if(!c)return null;var b=[];y(c,"array")&&y(c[0],"array")&&(b=a.path.clone(c));b.length||J(c).replace(ma,function(c,a,m){var e=[];a.toLowerCase();m.replace(Z,function(c,a){a&&e.push(+a)});b.push([a].concat(e))});b.toString=a.path.toString;return b};a._.svgTransform2string=d;a._.rgTransform=RegExp("^[a-z][\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*-?\\.?\\d",
"i");a._.transform2matrix=f;a._unit2px=b;a._.getSomeDefs=u;a._.getSomeSVG=p;a.select=function(c){return x(G.doc.querySelector(c))};a.selectAll=function(c){c=G.doc.querySelectorAll(c);for(var b=(a.set||Array)(),m=0;m<c.length;m++)b.push(x(c[m]));return b};setInterval(function(){for(var c in E)if(E[h](c)){var a=E[c],b=a.node;("svg"!=a.type&&!b.ownerSVGElement||"svg"==a.type&&(!b.parentNode||"ownerSVGElement"in b.parentNode&&!b.ownerSVGElement))&&delete E[c]}},1E4);(function(c){function m(c){function a(c,
b){var m=v(c.node,b);(m=(m=m&&m.match(d))&&m[2])&&"#"==m.charAt()&&(m=m.substring(1))&&(f[m]=(f[m]||[]).concat(function(a){var m={};m[b]=ca(a);v(c.node,m)}))}function b(c){var a=v(c.node,"xlink:href");a&&"#"==a.charAt()&&(a=a.substring(1))&&(f[a]=(f[a]||[]).concat(function(a){c.attr("xlink:href","#"+a)}))}var e=c.selectAll("*"),h,d=/^\s*url\(("|'|)(.*)\1\)\s*$/;c=[];for(var f={},l=0,E=e.length;l<E;l++){h=e[l];a(h,"fill");a(h,"stroke");a(h,"filter");a(h,"mask");a(h,"clip-path");b(h);var t=v(h.node,
"id");t&&(v(h.node,{id:h.id}),c.push({old:t,id:h.id}))}l=0;for(E=c.length;l<E;l++)if(e=f[c[l].old])for(h=0,t=e.length;h<t;h++)e[h](c[l].id)}function e(c,a,b){return function(m){m=m.slice(c,a);1==m.length&&(m=m[0]);return b?b(m):m}}function d(c){return function(){var a=c?"<"+this.type:"",b=this.node.attributes,m=this.node.childNodes;if(c)for(var e=0,h=b.length;e<h;e++)a+=" "+b[e].name+'="'+b[e].value.replace(/"/g,'\\"')+'"';if(m.length){c&&(a+=">");e=0;for(h=m.length;e<h;e++)3==m[e].nodeType?a+=m[e].nodeValue:
1==m[e].nodeType&&(a+=x(m[e]).toString());c&&(a+="</"+this.type+">")}else c&&(a+="/>");return a}}c.attr=function(c,a){if(!c)return this;if(y(c,"string"))if(1<arguments.length){var b={};b[c]=a;c=b}else return k("snap.util.getattr."+c,this).firstDefined();for(var m in c)c[h](m)&&k("snap.util.attr."+m,this,c[m]);return this};c.getBBox=function(c){if(!a.Matrix||!a.path)return this.node.getBBox();var b=this,m=new a.Matrix;if(b.removed)return a._.box();for(;"use"==b.type;)if(c||(m=m.add(b.transform().localMatrix.translate(b.attr("x")||
0,b.attr("y")||0))),b.original)b=b.original;else var e=b.attr("xlink:href"),b=b.original=b.node.ownerDocument.getElementById(e.substring(e.indexOf("#")+1));var e=b._,h=a.path.get[b.type]||a.path.get.deflt;try{if(c)return e.bboxwt=h?a.path.getBBox(b.realPath=h(b)):a._.box(b.node.getBBox()),a._.box(e.bboxwt);b.realPath=h(b);b.matrix=b.transform().localMatrix;e.bbox=a.path.getBBox(a.path.map(b.realPath,m.add(b.matrix)));return a._.box(e.bbox)}catch(d){return a._.box()}};var f=function(){return this.string};
c.transform=function(c){var b=this._;if(null==c){var m=this;c=new a.Matrix(this.node.getCTM());for(var e=n(this),h=[e],d=new a.Matrix,l=e.toTransformString(),b=J(e)==J(this.matrix)?J(b.transform):l;"svg"!=m.type&&(m=m.parent());)h.push(n(m));for(m=h.length;m--;)d.add(h[m]);return{string:b,globalMatrix:c,totalMatrix:d,localMatrix:e,diffMatrix:c.clone().add(e.invert()),global:c.toTransformString(),total:d.toTransformString(),local:l,toString:f}}c instanceof a.Matrix?this.matrix=c:n(this,c);this.node&&
("linearGradient"==this.type||"radialGradient"==this.type?v(this.node,{gradientTransform:this.matrix}):"pattern"==this.type?v(this.node,{patternTransform:this.matrix}):v(this.node,{transform:this.matrix}));return this};c.parent=function(){return x(this.node.parentNode)};c.append=c.add=function(c){if(c){if("set"==c.type){var a=this;c.forEach(function(c){a.add(c)});return this}c=x(c);this.node.appendChild(c.node);c.paper=this.paper}return this};c.appendTo=function(c){c&&(c=x(c),c.append(this));return this};
c.prepend=function(c){if(c){if("set"==c.type){var a=this,b;c.forEach(function(c){b?b.after(c):a.prepend(c);b=c});return this}c=x(c);var m=c.parent();this.node.insertBefore(c.node,this.node.firstChild);this.add&&this.add();c.paper=this.paper;this.parent()&&this.parent().add();m&&m.add()}return this};c.prependTo=function(c){c=x(c);c.prepend(this);return this};c.before=function(c){if("set"==c.type){var a=this;c.forEach(function(c){var b=c.parent();a.node.parentNode.insertBefore(c.node,a.node);b&&b.add()});
this.parent().add();return this}c=x(c);var b=c.parent();this.node.parentNode.insertBefore(c.node,this.node);this.parent()&&this.parent().add();b&&b.add();c.paper=this.paper;return this};c.after=function(c){c=x(c);var a=c.parent();this.node.nextSibling?this.node.parentNode.insertBefore(c.node,this.node.nextSibling):this.node.parentNode.appendChild(c.node);this.parent()&&this.parent().add();a&&a.add();c.paper=this.paper;return this};c.insertBefore=function(c){c=x(c);var a=this.parent();c.node.parentNode.insertBefore(this.node,
c.node);this.paper=c.paper;a&&a.add();c.parent()&&c.parent().add();return this};c.insertAfter=function(c){c=x(c);var a=this.parent();c.node.parentNode.insertBefore(this.node,c.node.nextSibling);this.paper=c.paper;a&&a.add();c.parent()&&c.parent().add();return this};c.remove=function(){var c=this.parent();this.node.parentNode&&this.node.parentNode.removeChild(this.node);delete this.paper;this.removed=!0;c&&c.add();return this};c.select=function(c){return x(this.node.querySelector(c))};c.selectAll=
function(c){c=this.node.querySelectorAll(c);for(var b=(a.set||Array)(),m=0;m<c.length;m++)b.push(x(c[m]));return b};c.asPX=function(c,a){null==a&&(a=this.attr(c));return+b(this,c,a)};c.use=function(){var c,a=this.node.id;a||(a=this.id,v(this.node,{id:a}));c="linearGradient"==this.type||"radialGradient"==this.type||"pattern"==this.type?r(this.type,this.node.parentNode):r("use",this.node.parentNode);v(c.node,{"xlink:href":"#"+a});c.original=this;return c};var l=/\S+/g;c.addClass=function(c){var a=(c||
"").match(l)||[];c=this.node;var b=c.className.baseVal,m=b.match(l)||[],e,h,d;if(a.length){for(e=0;d=a[e++];)h=m.indexOf(d),~h||m.push(d);a=m.join(" ");b!=a&&(c.className.baseVal=a)}return this};c.removeClass=function(c){var a=(c||"").match(l)||[];c=this.node;var b=c.className.baseVal,m=b.match(l)||[],e,h;if(m.length){for(e=0;h=a[e++];)h=m.indexOf(h),~h&&m.splice(h,1);a=m.join(" ");b!=a&&(c.className.baseVal=a)}return this};c.hasClass=function(c){return!!~(this.node.className.baseVal.match(l)||[]).indexOf(c)};
c.toggleClass=function(c,a){if(null!=a)return a?this.addClass(c):this.removeClass(c);var b=(c||"").match(l)||[],m=this.node,e=m.className.baseVal,h=e.match(l)||[],d,f,E;for(d=0;E=b[d++];)f=h.indexOf(E),~f?h.splice(f,1):h.push(E);b=h.join(" ");e!=b&&(m.className.baseVal=b);return this};c.clone=function(){var c=x(this.node.cloneNode(!0));v(c.node,"id")&&v(c.node,{id:c.id});m(c);c.insertAfter(this);return c};c.toDefs=function(){u(this).appendChild(this.node);return this};c.pattern=c.toPattern=function(c,
a,b,m){var e=r("pattern",u(this));null==c&&(c=this.getBBox());y(c,"object")&&"x"in c&&(a=c.y,b=c.width,m=c.height,c=c.x);v(e.node,{x:c,y:a,width:b,height:m,patternUnits:"userSpaceOnUse",id:e.id,viewBox:[c,a,b,m].join(" ")});e.node.appendChild(this.node);return e};c.marker=function(c,a,b,m,e,h){var d=r("marker",u(this));null==c&&(c=this.getBBox());y(c,"object")&&"x"in c&&(a=c.y,b=c.width,m=c.height,e=c.refX||c.cx,h=c.refY||c.cy,c=c.x);v(d.node,{viewBox:[c,a,b,m].join(" "),markerWidth:b,markerHeight:m,
orient:"auto",refX:e||0,refY:h||0,id:d.id});d.node.appendChild(this.node);return d};var E=function(c,a,b,m){"function"!=typeof b||b.length||(m=b,b=L.linear);this.attr=c;this.dur=a;b&&(this.easing=b);m&&(this.callback=m)};a._.Animation=E;a.animation=function(c,a,b,m){return new E(c,a,b,m)};c.inAnim=function(){var c=[],a;for(a in this.anims)this.anims[h](a)&&function(a){c.push({anim:new E(a._attrs,a.dur,a.easing,a._callback),mina:a,curStatus:a.status(),status:function(c){return a.status(c)},stop:function(){a.stop()}})}(this.anims[a]);
return c};a.animate=function(c,a,b,m,e,h){"function"!=typeof e||e.length||(h=e,e=L.linear);var d=L.time();c=L(c,a,d,d+m,L.time,b,e);h&&k.once("mina.finish."+c.id,h);return c};c.stop=function(){for(var c=this.inAnim(),a=0,b=c.length;a<b;a++)c[a].stop();return this};c.animate=function(c,a,b,m){"function"!=typeof b||b.length||(m=b,b=L.linear);c instanceof E&&(m=c.callback,b=c.easing,a=b.dur,c=c.attr);var d=[],f=[],l={},t,ca,n,T=this,q;for(q in c)if(c[h](q)){T.equal?(n=T.equal(q,J(c[q])),t=n.from,ca=
n.to,n=n.f):(t=+T.attr(q),ca=+c[q]);var la=y(t,"array")?t.length:1;l[q]=e(d.length,d.length+la,n);d=d.concat(t);f=f.concat(ca)}t=L.time();var p=L(d,f,t,t+a,L.time,function(c){var a={},b;for(b in l)l[h](b)&&(a[b]=l[b](c));T.attr(a)},b);T.anims[p.id]=p;p._attrs=c;p._callback=m;k("snap.animcreated."+T.id,p);k.once("mina.finish."+p.id,function(){delete T.anims[p.id];m&&m.call(T)});k.once("mina.stop."+p.id,function(){delete T.anims[p.id]});return T};var T={};c.data=function(c,b){var m=T[this.id]=T[this.id]||
{};if(0==arguments.length)return k("snap.data.get."+this.id,this,m,null),m;if(1==arguments.length){if(a.is(c,"object")){for(var e in c)c[h](e)&&this.data(e,c[e]);return this}k("snap.data.get."+this.id,this,m[c],c);return m[c]}m[c]=b;k("snap.data.set."+this.id,this,b,c);return this};c.removeData=function(c){null==c?T[this.id]={}:T[this.id]&&delete T[this.id][c];return this};c.outerSVG=c.toString=d(1);c.innerSVG=d()})(e.prototype);a.parse=function(c){var a=G.doc.createDocumentFragment(),b=!0,m=G.doc.createElement("div");
c=J(c);c.match(/^\s*<\s*svg(?:\s|>)/)||(c="<svg>"+c+"</svg>",b=!1);m.innerHTML=c;if(c=m.getElementsByTagName("svg")[0])if(b)a=c;else for(;c.firstChild;)a.appendChild(c.firstChild);m.innerHTML=aa;return new l(a)};l.prototype.select=e.prototype.select;l.prototype.selectAll=e.prototype.selectAll;a.fragment=function(){for(var c=Array.prototype.slice.call(arguments,0),b=G.doc.createDocumentFragment(),m=0,e=c.length;m<e;m++){var h=c[m];h.node&&h.node.nodeType&&b.appendChild(h.node);h.nodeType&&b.appendChild(h);
"string"==typeof h&&b.appendChild(a.parse(h).node)}return new l(b)};a._.make=r;a._.wrap=x;s.prototype.el=function(c,a){var b=r(c,this.node);a&&b.attr(a);return b};k.on("snap.util.getattr",function(){var c=k.nt(),c=c.substring(c.lastIndexOf(".")+1),a=c.replace(/[A-Z]/g,function(c){return"-"+c.toLowerCase()});return pa[h](a)?this.node.ownerDocument.defaultView.getComputedStyle(this.node,null).getPropertyValue(a):v(this.node,c)});var pa={"alignment-baseline":0,"baseline-shift":0,clip:0,"clip-path":0,
"clip-rule":0,color:0,"color-interpolation":0,"color-interpolation-filters":0,"color-profile":0,"color-rendering":0,cursor:0,direction:0,display:0,"dominant-baseline":0,"enable-background":0,fill:0,"fill-opacity":0,"fill-rule":0,filter:0,"flood-color":0,"flood-opacity":0,font:0,"font-family":0,"font-size":0,"font-size-adjust":0,"font-stretch":0,"font-style":0,"font-variant":0,"font-weight":0,"glyph-orientation-horizontal":0,"glyph-orientation-vertical":0,"image-rendering":0,kerning:0,"letter-spacing":0,
"lighting-color":0,marker:0,"marker-end":0,"marker-mid":0,"marker-start":0,mask:0,opacity:0,overflow:0,"pointer-events":0,"shape-rendering":0,"stop-color":0,"stop-opacity":0,stroke:0,"stroke-dasharray":0,"stroke-dashoffset":0,"stroke-linecap":0,"stroke-linejoin":0,"stroke-miterlimit":0,"stroke-opacity":0,"stroke-width":0,"text-anchor":0,"text-decoration":0,"text-rendering":0,"unicode-bidi":0,visibility:0,"word-spacing":0,"writing-mode":0};k.on("snap.util.attr",function(c){var a=k.nt(),b={},a=a.substring(a.lastIndexOf(".")+
1);b[a]=c;var m=a.replace(/-(\w)/gi,function(c,a){return a.toUpperCase()}),a=a.replace(/[A-Z]/g,function(c){return"-"+c.toLowerCase()});pa[h](a)?this.node.style[m]=null==c?aa:c:v(this.node,b)});a.ajax=function(c,a,b,m){var e=new XMLHttpRequest,h=V();if(e){if(y(a,"function"))m=b,b=a,a=null;else if(y(a,"object")){var d=[],f;for(f in a)a.hasOwnProperty(f)&&d.push(encodeURIComponent(f)+"="+encodeURIComponent(a[f]));a=d.join("&")}e.open(a?"POST":"GET",c,!0);a&&(e.setRequestHeader("X-Requested-With","XMLHttpRequest"),
e.setRequestHeader("Content-type","application/x-www-form-urlencoded"));b&&(k.once("snap.ajax."+h+".0",b),k.once("snap.ajax."+h+".200",b),k.once("snap.ajax."+h+".304",b));e.onreadystatechange=function(){4==e.readyState&&k("snap.ajax."+h+"."+e.status,m,e)};if(4==e.readyState)return e;e.send(a);return e}};a.load=function(c,b,m){a.ajax(c,function(c){c=a.parse(c.responseText);m?b.call(m,c):b(c)})};a.getElementByPoint=function(c,a){var b,m,e=G.doc.elementFromPoint(c,a);if(G.win.opera&&"svg"==e.tagName){b=
e;m=b.getBoundingClientRect();b=b.ownerDocument;var h=b.body,d=b.documentElement;b=m.top+(g.win.pageYOffset||d.scrollTop||h.scrollTop)-(d.clientTop||h.clientTop||0);m=m.left+(g.win.pageXOffset||d.scrollLeft||h.scrollLeft)-(d.clientLeft||h.clientLeft||0);h=e.createSVGRect();h.x=c-m;h.y=a-b;h.width=h.height=1;b=e.getIntersectionList(h,null);b.length&&(e=b[b.length-1])}return e?x(e):null};a.plugin=function(c){c(a,e,s,G,l)};return G.win.Snap=a}();C.plugin(function(a,k,y,M,A){function w(a,d,f,b,q,e){null==
d&&"[object SVGMatrix]"==z.call(a)?(this.a=a.a,this.b=a.b,this.c=a.c,this.d=a.d,this.e=a.e,this.f=a.f):null!=a?(this.a=+a,this.b=+d,this.c=+f,this.d=+b,this.e=+q,this.f=+e):(this.a=1,this.c=this.b=0,this.d=1,this.f=this.e=0)}var z=Object.prototype.toString,d=String,f=Math;(function(n){function k(a){return a[0]*a[0]+a[1]*a[1]}function p(a){var d=f.sqrt(k(a));a[0]&&(a[0]/=d);a[1]&&(a[1]/=d)}n.add=function(a,d,e,f,n,p){var k=[[],[],[] ],u=[[this.a,this.c,this.e],[this.b,this.d,this.f],[0,0,1] ];d=[[a,
e,n],[d,f,p],[0,0,1] ];a&&a instanceof w&&(d=[[a.a,a.c,a.e],[a.b,a.d,a.f],[0,0,1] ]);for(a=0;3>a;a++)for(e=0;3>e;e++){for(f=n=0;3>f;f++)n+=u[a][f]*d[f][e];k[a][e]=n}this.a=k[0][0];this.b=k[1][0];this.c=k[0][1];this.d=k[1][1];this.e=k[0][2];this.f=k[1][2];return this};n.invert=function(){var a=this.a*this.d-this.b*this.c;return new w(this.d/a,-this.b/a,-this.c/a,this.a/a,(this.c*this.f-this.d*this.e)/a,(this.b*this.e-this.a*this.f)/a)};n.clone=function(){return new w(this.a,this.b,this.c,this.d,this.e,
this.f)};n.translate=function(a,d){return this.add(1,0,0,1,a,d)};n.scale=function(a,d,e,f){null==d&&(d=a);(e||f)&&this.add(1,0,0,1,e,f);this.add(a,0,0,d,0,0);(e||f)&&this.add(1,0,0,1,-e,-f);return this};n.rotate=function(b,d,e){b=a.rad(b);d=d||0;e=e||0;var l=+f.cos(b).toFixed(9);b=+f.sin(b).toFixed(9);this.add(l,b,-b,l,d,e);return this.add(1,0,0,1,-d,-e)};n.x=function(a,d){return a*this.a+d*this.c+this.e};n.y=function(a,d){return a*this.b+d*this.d+this.f};n.get=function(a){return+this[d.fromCharCode(97+
a)].toFixed(4)};n.toString=function(){return"matrix("+[this.get(0),this.get(1),this.get(2),this.get(3),this.get(4),this.get(5)].join()+")"};n.offset=function(){return[this.e.toFixed(4),this.f.toFixed(4)]};n.determinant=function(){return this.a*this.d-this.b*this.c};n.split=function(){var b={};b.dx=this.e;b.dy=this.f;var d=[[this.a,this.c],[this.b,this.d] ];b.scalex=f.sqrt(k(d[0]));p(d[0]);b.shear=d[0][0]*d[1][0]+d[0][1]*d[1][1];d[1]=[d[1][0]-d[0][0]*b.shear,d[1][1]-d[0][1]*b.shear];b.scaley=f.sqrt(k(d[1]));
p(d[1]);b.shear/=b.scaley;0>this.determinant()&&(b.scalex=-b.scalex);var e=-d[0][1],d=d[1][1];0>d?(b.rotate=a.deg(f.acos(d)),0>e&&(b.rotate=360-b.rotate)):b.rotate=a.deg(f.asin(e));b.isSimple=!+b.shear.toFixed(9)&&(b.scalex.toFixed(9)==b.scaley.toFixed(9)||!b.rotate);b.isSuperSimple=!+b.shear.toFixed(9)&&b.scalex.toFixed(9)==b.scaley.toFixed(9)&&!b.rotate;b.noRotation=!+b.shear.toFixed(9)&&!b.rotate;return b};n.toTransformString=function(a){a=a||this.split();if(+a.shear.toFixed(9))return"m"+[this.get(0),
this.get(1),this.get(2),this.get(3),this.get(4),this.get(5)];a.scalex=+a.scalex.toFixed(4);a.scaley=+a.scaley.toFixed(4);a.rotate=+a.rotate.toFixed(4);return(a.dx||a.dy?"t"+[+a.dx.toFixed(4),+a.dy.toFixed(4)]:"")+(1!=a.scalex||1!=a.scaley?"s"+[a.scalex,a.scaley,0,0]:"")+(a.rotate?"r"+[+a.rotate.toFixed(4),0,0]:"")}})(w.prototype);a.Matrix=w;a.matrix=function(a,d,f,b,k,e){return new w(a,d,f,b,k,e)}});C.plugin(function(a,v,y,M,A){function w(h){return function(d){k.stop();d instanceof A&&1==d.node.childNodes.length&&
("radialGradient"==d.node.firstChild.tagName||"linearGradient"==d.node.firstChild.tagName||"pattern"==d.node.firstChild.tagName)&&(d=d.node.firstChild,b(this).appendChild(d),d=u(d));if(d instanceof v)if("radialGradient"==d.type||"linearGradient"==d.type||"pattern"==d.type){d.node.id||e(d.node,{id:d.id});var f=l(d.node.id)}else f=d.attr(h);else f=a.color(d),f.error?(f=a(b(this).ownerSVGElement).gradient(d))?(f.node.id||e(f.node,{id:f.id}),f=l(f.node.id)):f=d:f=r(f);d={};d[h]=f;e(this.node,d);this.node.style[h]=
x}}function z(a){k.stop();a==+a&&(a+="px");this.node.style.fontSize=a}function d(a){var b=[];a=a.childNodes;for(var e=0,f=a.length;e<f;e++){var l=a[e];3==l.nodeType&&b.push(l.nodeValue);"tspan"==l.tagName&&(1==l.childNodes.length&&3==l.firstChild.nodeType?b.push(l.firstChild.nodeValue):b.push(d(l)))}return b}function f(){k.stop();return this.node.style.fontSize}var n=a._.make,u=a._.wrap,p=a.is,b=a._.getSomeDefs,q=/^url\(#?([^)]+)\)$/,e=a._.$,l=a.url,r=String,s=a._.separator,x="";k.on("snap.util.attr.mask",
function(a){if(a instanceof v||a instanceof A){k.stop();a instanceof A&&1==a.node.childNodes.length&&(a=a.node.firstChild,b(this).appendChild(a),a=u(a));if("mask"==a.type)var d=a;else d=n("mask",b(this)),d.node.appendChild(a.node);!d.node.id&&e(d.node,{id:d.id});e(this.node,{mask:l(d.id)})}});(function(a){k.on("snap.util.attr.clip",a);k.on("snap.util.attr.clip-path",a);k.on("snap.util.attr.clipPath",a)})(function(a){if(a instanceof v||a instanceof A){k.stop();if("clipPath"==a.type)var d=a;else d=
n("clipPath",b(this)),d.node.appendChild(a.node),!d.node.id&&e(d.node,{id:d.id});e(this.node,{"clip-path":l(d.id)})}});k.on("snap.util.attr.fill",w("fill"));k.on("snap.util.attr.stroke",w("stroke"));var G=/^([lr])(?:\(([^)]*)\))?(.*)$/i;k.on("snap.util.grad.parse",function(a){a=r(a);var b=a.match(G);if(!b)return null;a=b[1];var e=b[2],b=b[3],e=e.split(/\s*,\s*/).map(function(a){return+a==a?+a:a});1==e.length&&0==e[0]&&(e=[]);b=b.split("-");b=b.map(function(a){a=a.split(":");var b={color:a[0]};a[1]&&
(b.offset=parseFloat(a[1]));return b});return{type:a,params:e,stops:b}});k.on("snap.util.attr.d",function(b){k.stop();p(b,"array")&&p(b[0],"array")&&(b=a.path.toString.call(b));b=r(b);b.match(/[ruo]/i)&&(b=a.path.toAbsolute(b));e(this.node,{d:b})})(-1);k.on("snap.util.attr.#text",function(a){k.stop();a=r(a);for(a=M.doc.createTextNode(a);this.node.firstChild;)this.node.removeChild(this.node.firstChild);this.node.appendChild(a)})(-1);k.on("snap.util.attr.path",function(a){k.stop();this.attr({d:a})})(-1);
k.on("snap.util.attr.class",function(a){k.stop();this.node.className.baseVal=a})(-1);k.on("snap.util.attr.viewBox",function(a){a=p(a,"object")&&"x"in a?[a.x,a.y,a.width,a.height].join(" "):p(a,"array")?a.join(" "):a;e(this.node,{viewBox:a});k.stop()})(-1);k.on("snap.util.attr.transform",function(a){this.transform(a);k.stop()})(-1);k.on("snap.util.attr.r",function(a){"rect"==this.type&&(k.stop(),e(this.node,{rx:a,ry:a}))})(-1);k.on("snap.util.attr.textpath",function(a){k.stop();if("text"==this.type){var d,
f;if(!a&&this.textPath){for(a=this.textPath;a.node.firstChild;)this.node.appendChild(a.node.firstChild);a.remove();delete this.textPath}else if(p(a,"string")?(d=b(this),a=u(d.parentNode).path(a),d.appendChild(a.node),d=a.id,a.attr({id:d})):(a=u(a),a instanceof v&&(d=a.attr("id"),d||(d=a.id,a.attr({id:d})))),d)if(a=this.textPath,f=this.node,a)a.attr({"xlink:href":"#"+d});else{for(a=e("textPath",{"xlink:href":"#"+d});f.firstChild;)a.appendChild(f.firstChild);f.appendChild(a);this.textPath=u(a)}}})(-1);
k.on("snap.util.attr.text",function(a){if("text"==this.type){for(var b=this.node,d=function(a){var b=e("tspan");if(p(a,"array"))for(var f=0;f<a.length;f++)b.appendChild(d(a[f]));else b.appendChild(M.doc.createTextNode(a));b.normalize&&b.normalize();return b};b.firstChild;)b.removeChild(b.firstChild);for(a=d(a);a.firstChild;)b.appendChild(a.firstChild)}k.stop()})(-1);k.on("snap.util.attr.fontSize",z)(-1);k.on("snap.util.attr.font-size",z)(-1);k.on("snap.util.getattr.transform",function(){k.stop();
return this.transform()})(-1);k.on("snap.util.getattr.textpath",function(){k.stop();return this.textPath})(-1);(function(){function b(d){return function(){k.stop();var b=M.doc.defaultView.getComputedStyle(this.node,null).getPropertyValue("marker-"+d);return"none"==b?b:a(M.doc.getElementById(b.match(q)[1]))}}function d(a){return function(b){k.stop();var d="marker"+a.charAt(0).toUpperCase()+a.substring(1);if(""==b||!b)this.node.style[d]="none";else if("marker"==b.type){var f=b.node.id;f||e(b.node,{id:b.id});
this.node.style[d]=l(f)}}}k.on("snap.util.getattr.marker-end",b("end"))(-1);k.on("snap.util.getattr.markerEnd",b("end"))(-1);k.on("snap.util.getattr.marker-start",b("start"))(-1);k.on("snap.util.getattr.markerStart",b("start"))(-1);k.on("snap.util.getattr.marker-mid",b("mid"))(-1);k.on("snap.util.getattr.markerMid",b("mid"))(-1);k.on("snap.util.attr.marker-end",d("end"))(-1);k.on("snap.util.attr.markerEnd",d("end"))(-1);k.on("snap.util.attr.marker-start",d("start"))(-1);k.on("snap.util.attr.markerStart",
d("start"))(-1);k.on("snap.util.attr.marker-mid",d("mid"))(-1);k.on("snap.util.attr.markerMid",d("mid"))(-1)})();k.on("snap.util.getattr.r",function(){if("rect"==this.type&&e(this.node,"rx")==e(this.node,"ry"))return k.stop(),e(this.node,"rx")})(-1);k.on("snap.util.getattr.text",function(){if("text"==this.type||"tspan"==this.type){k.stop();var a=d(this.node);return 1==a.length?a[0]:a}})(-1);k.on("snap.util.getattr.#text",function(){return this.node.textContent})(-1);k.on("snap.util.getattr.viewBox",
function(){k.stop();var b=e(this.node,"viewBox");if(b)return b=b.split(s),a._.box(+b[0],+b[1],+b[2],+b[3])})(-1);k.on("snap.util.getattr.points",function(){var a=e(this.node,"points");k.stop();if(a)return a.split(s)})(-1);k.on("snap.util.getattr.path",function(){var a=e(this.node,"d");k.stop();return a})(-1);k.on("snap.util.getattr.class",function(){return this.node.className.baseVal})(-1);k.on("snap.util.getattr.fontSize",f)(-1);k.on("snap.util.getattr.font-size",f)(-1)});C.plugin(function(a,v,y,
M,A){function w(a){return a}function z(a){return function(b){return+b.toFixed(3)+a}}var d={"+":function(a,b){return a+b},"-":function(a,b){return a-b},"/":function(a,b){return a/b},"*":function(a,b){return a*b}},f=String,n=/[a-z]+$/i,u=/^\s*([+\-\/*])\s*=\s*([\d.eE+\-]+)\s*([^\d\s]+)?\s*$/;k.on("snap.util.attr",function(a){if(a=f(a).match(u)){var b=k.nt(),b=b.substring(b.lastIndexOf(".")+1),q=this.attr(b),e={};k.stop();var l=a[3]||"",r=q.match(n),s=d[a[1] ];r&&r==l?a=s(parseFloat(q),+a[2]):(q=this.asPX(b),
a=s(this.asPX(b),this.asPX(b,a[2]+l)));isNaN(q)||isNaN(a)||(e[b]=a,this.attr(e))}})(-10);k.on("snap.util.equal",function(a,b){var q=f(this.attr(a)||""),e=f(b).match(u);if(e){k.stop();var l=e[3]||"",r=q.match(n),s=d[e[1] ];if(r&&r==l)return{from:parseFloat(q),to:s(parseFloat(q),+e[2]),f:z(r)};q=this.asPX(a);return{from:q,to:s(q,this.asPX(a,e[2]+l)),f:w}}})(-10)});C.plugin(function(a,v,y,M,A){var w=y.prototype,z=a.is;w.rect=function(a,d,k,p,b,q){var e;null==q&&(q=b);z(a,"object")&&"[object Object]"==
a?e=a:null!=a&&(e={x:a,y:d,width:k,height:p},null!=b&&(e.rx=b,e.ry=q));return this.el("rect",e)};w.circle=function(a,d,k){var p;z(a,"object")&&"[object Object]"==a?p=a:null!=a&&(p={cx:a,cy:d,r:k});return this.el("circle",p)};var d=function(){function a(){this.parentNode.removeChild(this)}return function(d,k){var p=M.doc.createElement("img"),b=M.doc.body;p.style.cssText="position:absolute;left:-9999em;top:-9999em";p.onload=function(){k.call(p);p.onload=p.onerror=null;b.removeChild(p)};p.onerror=a;
b.appendChild(p);p.src=d}}();w.image=function(f,n,k,p,b){var q=this.el("image");if(z(f,"object")&&"src"in f)q.attr(f);else if(null!=f){var e={"xlink:href":f,preserveAspectRatio:"none"};null!=n&&null!=k&&(e.x=n,e.y=k);null!=p&&null!=b?(e.width=p,e.height=b):d(f,function(){a._.$(q.node,{width:this.offsetWidth,height:this.offsetHeight})});a._.$(q.node,e)}return q};w.ellipse=function(a,d,k,p){var b;z(a,"object")&&"[object Object]"==a?b=a:null!=a&&(b={cx:a,cy:d,rx:k,ry:p});return this.el("ellipse",b)};
w.path=function(a){var d;z(a,"object")&&!z(a,"array")?d=a:a&&(d={d:a});return this.el("path",d)};w.group=w.g=function(a){var d=this.el("g");1==arguments.length&&a&&!a.type?d.attr(a):arguments.length&&d.add(Array.prototype.slice.call(arguments,0));return d};w.svg=function(a,d,k,p,b,q,e,l){var r={};z(a,"object")&&null==d?r=a:(null!=a&&(r.x=a),null!=d&&(r.y=d),null!=k&&(r.width=k),null!=p&&(r.height=p),null!=b&&null!=q&&null!=e&&null!=l&&(r.viewBox=[b,q,e,l]));return this.el("svg",r)};w.mask=function(a){var d=
this.el("mask");1==arguments.length&&a&&!a.type?d.attr(a):arguments.length&&d.add(Array.prototype.slice.call(arguments,0));return d};w.ptrn=function(a,d,k,p,b,q,e,l){if(z(a,"object"))var r=a;else arguments.length?(r={},null!=a&&(r.x=a),null!=d&&(r.y=d),null!=k&&(r.width=k),null!=p&&(r.height=p),null!=b&&null!=q&&null!=e&&null!=l&&(r.viewBox=[b,q,e,l])):r={patternUnits:"userSpaceOnUse"};return this.el("pattern",r)};w.use=function(a){return null!=a?(make("use",this.node),a instanceof v&&(a.attr("id")||
a.attr({id:ID()}),a=a.attr("id")),this.el("use",{"xlink:href":a})):v.prototype.use.call(this)};w.text=function(a,d,k){var p={};z(a,"object")?p=a:null!=a&&(p={x:a,y:d,text:k||""});return this.el("text",p)};w.line=function(a,d,k,p){var b={};z(a,"object")?b=a:null!=a&&(b={x1:a,x2:k,y1:d,y2:p});return this.el("line",b)};w.polyline=function(a){1<arguments.length&&(a=Array.prototype.slice.call(arguments,0));var d={};z(a,"object")&&!z(a,"array")?d=a:null!=a&&(d={points:a});return this.el("polyline",d)};
w.polygon=function(a){1<arguments.length&&(a=Array.prototype.slice.call(arguments,0));var d={};z(a,"object")&&!z(a,"array")?d=a:null!=a&&(d={points:a});return this.el("polygon",d)};(function(){function d(){return this.selectAll("stop")}function n(b,d){var f=e("stop"),k={offset:+d+"%"};b=a.color(b);k["stop-color"]=b.hex;1>b.opacity&&(k["stop-opacity"]=b.opacity);e(f,k);this.node.appendChild(f);return this}function u(){if("linearGradient"==this.type){var b=e(this.node,"x1")||0,d=e(this.node,"x2")||
1,f=e(this.node,"y1")||0,k=e(this.node,"y2")||0;return a._.box(b,f,math.abs(d-b),math.abs(k-f))}b=this.node.r||0;return a._.box((this.node.cx||0.5)-b,(this.node.cy||0.5)-b,2*b,2*b)}function p(a,d){function f(a,b){for(var d=(b-u)/(a-w),e=w;e<a;e++)h[e].offset=+(+u+d*(e-w)).toFixed(2);w=a;u=b}var n=k("snap.util.grad.parse",null,d).firstDefined(),p;if(!n)return null;n.params.unshift(a);p="l"==n.type.toLowerCase()?b.apply(0,n.params):q.apply(0,n.params);n.type!=n.type.toLowerCase()&&e(p.node,{gradientUnits:"userSpaceOnUse"});
var h=n.stops,n=h.length,u=0,w=0;n--;for(var v=0;v<n;v++)"offset"in h[v]&&f(v,h[v].offset);h[n].offset=h[n].offset||100;f(n,h[n].offset);for(v=0;v<=n;v++){var y=h[v];p.addStop(y.color,y.offset)}return p}function b(b,k,p,q,w){b=a._.make("linearGradient",b);b.stops=d;b.addStop=n;b.getBBox=u;null!=k&&e(b.node,{x1:k,y1:p,x2:q,y2:w});return b}function q(b,k,p,q,w,h){b=a._.make("radialGradient",b);b.stops=d;b.addStop=n;b.getBBox=u;null!=k&&e(b.node,{cx:k,cy:p,r:q});null!=w&&null!=h&&e(b.node,{fx:w,fy:h});
return b}var e=a._.$;w.gradient=function(a){return p(this.defs,a)};w.gradientLinear=function(a,d,e,f){return b(this.defs,a,d,e,f)};w.gradientRadial=function(a,b,d,e,f){return q(this.defs,a,b,d,e,f)};w.toString=function(){var b=this.node.ownerDocument,d=b.createDocumentFragment(),b=b.createElement("div"),e=this.node.cloneNode(!0);d.appendChild(b);b.appendChild(e);a._.$(e,{xmlns:"http://www.w3.org/2000/svg"});b=b.innerHTML;d.removeChild(d.firstChild);return b};w.clear=function(){for(var a=this.node.firstChild,
b;a;)b=a.nextSibling,"defs"!=a.tagName?a.parentNode.removeChild(a):w.clear.call({node:a}),a=b}})()});C.plugin(function(a,k,y,M){function A(a){var b=A.ps=A.ps||{};b[a]?b[a].sleep=100:b[a]={sleep:100};setTimeout(function(){for(var d in b)b[L](d)&&d!=a&&(b[d].sleep--,!b[d].sleep&&delete b[d])});return b[a]}function w(a,b,d,e){null==a&&(a=b=d=e=0);null==b&&(b=a.y,d=a.width,e=a.height,a=a.x);return{x:a,y:b,width:d,w:d,height:e,h:e,x2:a+d,y2:b+e,cx:a+d/2,cy:b+e/2,r1:F.min(d,e)/2,r2:F.max(d,e)/2,r0:F.sqrt(d*
d+e*e)/2,path:s(a,b,d,e),vb:[a,b,d,e].join(" ")}}function z(){return this.join(",").replace(N,"$1")}function d(a){a=C(a);a.toString=z;return a}function f(a,b,d,h,f,k,l,n,p){if(null==p)return e(a,b,d,h,f,k,l,n);if(0>p||e(a,b,d,h,f,k,l,n)<p)p=void 0;else{var q=0.5,O=1-q,s;for(s=e(a,b,d,h,f,k,l,n,O);0.01<Z(s-p);)q/=2,O+=(s<p?1:-1)*q,s=e(a,b,d,h,f,k,l,n,O);p=O}return u(a,b,d,h,f,k,l,n,p)}function n(b,d){function e(a){return+(+a).toFixed(3)}return a._.cacher(function(a,h,l){a instanceof k&&(a=a.attr("d"));
a=I(a);for(var n,p,D,q,O="",s={},c=0,t=0,r=a.length;t<r;t++){D=a[t];if("M"==D[0])n=+D[1],p=+D[2];else{q=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6]);if(c+q>h){if(d&&!s.start){n=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6],h-c);O+=["C"+e(n.start.x),e(n.start.y),e(n.m.x),e(n.m.y),e(n.x),e(n.y)];if(l)return O;s.start=O;O=["M"+e(n.x),e(n.y)+"C"+e(n.n.x),e(n.n.y),e(n.end.x),e(n.end.y),e(D[5]),e(D[6])].join();c+=q;n=+D[5];p=+D[6];continue}if(!b&&!d)return n=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6],h-c)}c+=q;n=+D[5];p=+D[6]}O+=
D.shift()+D}s.end=O;return n=b?c:d?s:u(n,p,D[0],D[1],D[2],D[3],D[4],D[5],1)},null,a._.clone)}function u(a,b,d,e,h,f,k,l,n){var p=1-n,q=ma(p,3),s=ma(p,2),c=n*n,t=c*n,r=q*a+3*s*n*d+3*p*n*n*h+t*k,q=q*b+3*s*n*e+3*p*n*n*f+t*l,s=a+2*n*(d-a)+c*(h-2*d+a),t=b+2*n*(e-b)+c*(f-2*e+b),x=d+2*n*(h-d)+c*(k-2*h+d),c=e+2*n*(f-e)+c*(l-2*f+e);a=p*a+n*d;b=p*b+n*e;h=p*h+n*k;f=p*f+n*l;l=90-180*F.atan2(s-x,t-c)/S;return{x:r,y:q,m:{x:s,y:t},n:{x:x,y:c},start:{x:a,y:b},end:{x:h,y:f},alpha:l}}function p(b,d,e,h,f,n,k,l){a.is(b,
"array")||(b=[b,d,e,h,f,n,k,l]);b=U.apply(null,b);return w(b.min.x,b.min.y,b.max.x-b.min.x,b.max.y-b.min.y)}function b(a,b,d){return b>=a.x&&b<=a.x+a.width&&d>=a.y&&d<=a.y+a.height}function q(a,d){a=w(a);d=w(d);return b(d,a.x,a.y)||b(d,a.x2,a.y)||b(d,a.x,a.y2)||b(d,a.x2,a.y2)||b(a,d.x,d.y)||b(a,d.x2,d.y)||b(a,d.x,d.y2)||b(a,d.x2,d.y2)||(a.x<d.x2&&a.x>d.x||d.x<a.x2&&d.x>a.x)&&(a.y<d.y2&&a.y>d.y||d.y<a.y2&&d.y>a.y)}function e(a,b,d,e,h,f,n,k,l){null==l&&(l=1);l=(1<l?1:0>l?0:l)/2;for(var p=[-0.1252,
0.1252,-0.3678,0.3678,-0.5873,0.5873,-0.7699,0.7699,-0.9041,0.9041,-0.9816,0.9816],q=[0.2491,0.2491,0.2335,0.2335,0.2032,0.2032,0.1601,0.1601,0.1069,0.1069,0.0472,0.0472],s=0,c=0;12>c;c++)var t=l*p[c]+l,r=t*(t*(-3*a+9*d-9*h+3*n)+6*a-12*d+6*h)-3*a+3*d,t=t*(t*(-3*b+9*e-9*f+3*k)+6*b-12*e+6*f)-3*b+3*e,s=s+q[c]*F.sqrt(r*r+t*t);return l*s}function l(a,b,d){a=I(a);b=I(b);for(var h,f,l,n,k,s,r,O,x,c,t=d?0:[],w=0,v=a.length;w<v;w++)if(x=a[w],"M"==x[0])h=k=x[1],f=s=x[2];else{"C"==x[0]?(x=[h,f].concat(x.slice(1)),
h=x[6],f=x[7]):(x=[h,f,h,f,k,s,k,s],h=k,f=s);for(var G=0,y=b.length;G<y;G++)if(c=b[G],"M"==c[0])l=r=c[1],n=O=c[2];else{"C"==c[0]?(c=[l,n].concat(c.slice(1)),l=c[6],n=c[7]):(c=[l,n,l,n,r,O,r,O],l=r,n=O);var z;var K=x,B=c;z=d;var H=p(K),J=p(B);if(q(H,J)){for(var H=e.apply(0,K),J=e.apply(0,B),H=~~(H/8),J=~~(J/8),U=[],A=[],F={},M=z?0:[],P=0;P<H+1;P++){var C=u.apply(0,K.concat(P/H));U.push({x:C.x,y:C.y,t:P/H})}for(P=0;P<J+1;P++)C=u.apply(0,B.concat(P/J)),A.push({x:C.x,y:C.y,t:P/J});for(P=0;P<H;P++)for(K=
0;K<J;K++){var Q=U[P],L=U[P+1],B=A[K],C=A[K+1],N=0.001>Z(L.x-Q.x)?"y":"x",S=0.001>Z(C.x-B.x)?"y":"x",R;R=Q.x;var Y=Q.y,V=L.x,ea=L.y,fa=B.x,ga=B.y,ha=C.x,ia=C.y;if(W(R,V)<X(fa,ha)||X(R,V)>W(fa,ha)||W(Y,ea)<X(ga,ia)||X(Y,ea)>W(ga,ia))R=void 0;else{var $=(R*ea-Y*V)*(fa-ha)-(R-V)*(fa*ia-ga*ha),aa=(R*ea-Y*V)*(ga-ia)-(Y-ea)*(fa*ia-ga*ha),ja=(R-V)*(ga-ia)-(Y-ea)*(fa-ha);if(ja){var $=$/ja,aa=aa/ja,ja=+$.toFixed(2),ba=+aa.toFixed(2);R=ja<+X(R,V).toFixed(2)||ja>+W(R,V).toFixed(2)||ja<+X(fa,ha).toFixed(2)||
ja>+W(fa,ha).toFixed(2)||ba<+X(Y,ea).toFixed(2)||ba>+W(Y,ea).toFixed(2)||ba<+X(ga,ia).toFixed(2)||ba>+W(ga,ia).toFixed(2)?void 0:{x:$,y:aa}}else R=void 0}R&&F[R.x.toFixed(4)]!=R.y.toFixed(4)&&(F[R.x.toFixed(4)]=R.y.toFixed(4),Q=Q.t+Z((R[N]-Q[N])/(L[N]-Q[N]))*(L.t-Q.t),B=B.t+Z((R[S]-B[S])/(C[S]-B[S]))*(C.t-B.t),0<=Q&&1>=Q&&0<=B&&1>=B&&(z?M++:M.push({x:R.x,y:R.y,t1:Q,t2:B})))}z=M}else z=z?0:[];if(d)t+=z;else{H=0;for(J=z.length;H<J;H++)z[H].segment1=w,z[H].segment2=G,z[H].bez1=x,z[H].bez2=c;t=t.concat(z)}}}return t}
function r(a){var b=A(a);if(b.bbox)return C(b.bbox);if(!a)return w();a=I(a);for(var d=0,e=0,h=[],f=[],l,n=0,k=a.length;n<k;n++)l=a[n],"M"==l[0]?(d=l[1],e=l[2],h.push(d),f.push(e)):(d=U(d,e,l[1],l[2],l[3],l[4],l[5],l[6]),h=h.concat(d.min.x,d.max.x),f=f.concat(d.min.y,d.max.y),d=l[5],e=l[6]);a=X.apply(0,h);l=X.apply(0,f);h=W.apply(0,h);f=W.apply(0,f);f=w(a,l,h-a,f-l);b.bbox=C(f);return f}function s(a,b,d,e,h){if(h)return[["M",+a+ +h,b],["l",d-2*h,0],["a",h,h,0,0,1,h,h],["l",0,e-2*h],["a",h,h,0,0,1,
-h,h],["l",2*h-d,0],["a",h,h,0,0,1,-h,-h],["l",0,2*h-e],["a",h,h,0,0,1,h,-h],["z"] ];a=[["M",a,b],["l",d,0],["l",0,e],["l",-d,0],["z"] ];a.toString=z;return a}function x(a,b,d,e,h){null==h&&null==e&&(e=d);a=+a;b=+b;d=+d;e=+e;if(null!=h){var f=Math.PI/180,l=a+d*Math.cos(-e*f);a+=d*Math.cos(-h*f);var n=b+d*Math.sin(-e*f);b+=d*Math.sin(-h*f);d=[["M",l,n],["A",d,d,0,+(180<h-e),0,a,b] ]}else d=[["M",a,b],["m",0,-e],["a",d,e,0,1,1,0,2*e],["a",d,e,0,1,1,0,-2*e],["z"] ];d.toString=z;return d}function G(b){var e=
A(b);if(e.abs)return d(e.abs);Q(b,"array")&&Q(b&&b[0],"array")||(b=a.parsePathString(b));if(!b||!b.length)return[["M",0,0] ];var h=[],f=0,l=0,n=0,k=0,p=0;"M"==b[0][0]&&(f=+b[0][1],l=+b[0][2],n=f,k=l,p++,h[0]=["M",f,l]);for(var q=3==b.length&&"M"==b[0][0]&&"R"==b[1][0].toUpperCase()&&"Z"==b[2][0].toUpperCase(),s,r,w=p,c=b.length;w<c;w++){h.push(s=[]);r=b[w];p=r[0];if(p!=p.toUpperCase())switch(s[0]=p.toUpperCase(),s[0]){case "A":s[1]=r[1];s[2]=r[2];s[3]=r[3];s[4]=r[4];s[5]=r[5];s[6]=+r[6]+f;s[7]=+r[7]+
l;break;case "V":s[1]=+r[1]+l;break;case "H":s[1]=+r[1]+f;break;case "R":for(var t=[f,l].concat(r.slice(1)),u=2,v=t.length;u<v;u++)t[u]=+t[u]+f,t[++u]=+t[u]+l;h.pop();h=h.concat(P(t,q));break;case "O":h.pop();t=x(f,l,r[1],r[2]);t.push(t[0]);h=h.concat(t);break;case "U":h.pop();h=h.concat(x(f,l,r[1],r[2],r[3]));s=["U"].concat(h[h.length-1].slice(-2));break;case "M":n=+r[1]+f,k=+r[2]+l;default:for(u=1,v=r.length;u<v;u++)s[u]=+r[u]+(u%2?f:l)}else if("R"==p)t=[f,l].concat(r.slice(1)),h.pop(),h=h.concat(P(t,
q)),s=["R"].concat(r.slice(-2));else if("O"==p)h.pop(),t=x(f,l,r[1],r[2]),t.push(t[0]),h=h.concat(t);else if("U"==p)h.pop(),h=h.concat(x(f,l,r[1],r[2],r[3])),s=["U"].concat(h[h.length-1].slice(-2));else for(t=0,u=r.length;t<u;t++)s[t]=r[t];p=p.toUpperCase();if("O"!=p)switch(s[0]){case "Z":f=+n;l=+k;break;case "H":f=s[1];break;case "V":l=s[1];break;case "M":n=s[s.length-2],k=s[s.length-1];default:f=s[s.length-2],l=s[s.length-1]}}h.toString=z;e.abs=d(h);return h}function h(a,b,d,e){return[a,b,d,e,d,
e]}function J(a,b,d,e,h,f){var l=1/3,n=2/3;return[l*a+n*d,l*b+n*e,l*h+n*d,l*f+n*e,h,f]}function K(b,d,e,h,f,l,n,k,p,s){var r=120*S/180,q=S/180*(+f||0),c=[],t,x=a._.cacher(function(a,b,c){var d=a*F.cos(c)-b*F.sin(c);a=a*F.sin(c)+b*F.cos(c);return{x:d,y:a}});if(s)v=s[0],t=s[1],l=s[2],u=s[3];else{t=x(b,d,-q);b=t.x;d=t.y;t=x(k,p,-q);k=t.x;p=t.y;F.cos(S/180*f);F.sin(S/180*f);t=(b-k)/2;v=(d-p)/2;u=t*t/(e*e)+v*v/(h*h);1<u&&(u=F.sqrt(u),e*=u,h*=u);var u=e*e,w=h*h,u=(l==n?-1:1)*F.sqrt(Z((u*w-u*v*v-w*t*t)/
(u*v*v+w*t*t)));l=u*e*v/h+(b+k)/2;var u=u*-h*t/e+(d+p)/2,v=F.asin(((d-u)/h).toFixed(9));t=F.asin(((p-u)/h).toFixed(9));v=b<l?S-v:v;t=k<l?S-t:t;0>v&&(v=2*S+v);0>t&&(t=2*S+t);n&&v>t&&(v-=2*S);!n&&t>v&&(t-=2*S)}if(Z(t-v)>r){var c=t,w=k,G=p;t=v+r*(n&&t>v?1:-1);k=l+e*F.cos(t);p=u+h*F.sin(t);c=K(k,p,e,h,f,0,n,w,G,[t,c,l,u])}l=t-v;f=F.cos(v);r=F.sin(v);n=F.cos(t);t=F.sin(t);l=F.tan(l/4);e=4/3*e*l;l*=4/3*h;h=[b,d];b=[b+e*r,d-l*f];d=[k+e*t,p-l*n];k=[k,p];b[0]=2*h[0]-b[0];b[1]=2*h[1]-b[1];if(s)return[b,d,k].concat(c);
c=[b,d,k].concat(c).join().split(",");s=[];k=0;for(p=c.length;k<p;k++)s[k]=k%2?x(c[k-1],c[k],q).y:x(c[k],c[k+1],q).x;return s}function U(a,b,d,e,h,f,l,k){for(var n=[],p=[[],[] ],s,r,c,t,q=0;2>q;++q)0==q?(r=6*a-12*d+6*h,s=-3*a+9*d-9*h+3*l,c=3*d-3*a):(r=6*b-12*e+6*f,s=-3*b+9*e-9*f+3*k,c=3*e-3*b),1E-12>Z(s)?1E-12>Z(r)||(s=-c/r,0<s&&1>s&&n.push(s)):(t=r*r-4*c*s,c=F.sqrt(t),0>t||(t=(-r+c)/(2*s),0<t&&1>t&&n.push(t),s=(-r-c)/(2*s),0<s&&1>s&&n.push(s)));for(r=q=n.length;q--;)s=n[q],c=1-s,p[0][q]=c*c*c*a+3*
c*c*s*d+3*c*s*s*h+s*s*s*l,p[1][q]=c*c*c*b+3*c*c*s*e+3*c*s*s*f+s*s*s*k;p[0][r]=a;p[1][r]=b;p[0][r+1]=l;p[1][r+1]=k;p[0].length=p[1].length=r+2;return{min:{x:X.apply(0,p[0]),y:X.apply(0,p[1])},max:{x:W.apply(0,p[0]),y:W.apply(0,p[1])}}}function I(a,b){var e=!b&&A(a);if(!b&&e.curve)return d(e.curve);var f=G(a),l=b&&G(b),n={x:0,y:0,bx:0,by:0,X:0,Y:0,qx:null,qy:null},k={x:0,y:0,bx:0,by:0,X:0,Y:0,qx:null,qy:null},p=function(a,b,c){if(!a)return["C",b.x,b.y,b.x,b.y,b.x,b.y];a[0]in{T:1,Q:1}||(b.qx=b.qy=null);
switch(a[0]){case "M":b.X=a[1];b.Y=a[2];break;case "A":a=["C"].concat(K.apply(0,[b.x,b.y].concat(a.slice(1))));break;case "S":"C"==c||"S"==c?(c=2*b.x-b.bx,b=2*b.y-b.by):(c=b.x,b=b.y);a=["C",c,b].concat(a.slice(1));break;case "T":"Q"==c||"T"==c?(b.qx=2*b.x-b.qx,b.qy=2*b.y-b.qy):(b.qx=b.x,b.qy=b.y);a=["C"].concat(J(b.x,b.y,b.qx,b.qy,a[1],a[2]));break;case "Q":b.qx=a[1];b.qy=a[2];a=["C"].concat(J(b.x,b.y,a[1],a[2],a[3],a[4]));break;case "L":a=["C"].concat(h(b.x,b.y,a[1],a[2]));break;case "H":a=["C"].concat(h(b.x,
b.y,a[1],b.y));break;case "V":a=["C"].concat(h(b.x,b.y,b.x,a[1]));break;case "Z":a=["C"].concat(h(b.x,b.y,b.X,b.Y))}return a},s=function(a,b){if(7<a[b].length){a[b].shift();for(var c=a[b];c.length;)q[b]="A",l&&(u[b]="A"),a.splice(b++,0,["C"].concat(c.splice(0,6)));a.splice(b,1);v=W(f.length,l&&l.length||0)}},r=function(a,b,c,d,e){a&&b&&"M"==a[e][0]&&"M"!=b[e][0]&&(b.splice(e,0,["M",d.x,d.y]),c.bx=0,c.by=0,c.x=a[e][1],c.y=a[e][2],v=W(f.length,l&&l.length||0))},q=[],u=[],c="",t="",x=0,v=W(f.length,
l&&l.length||0);for(;x<v;x++){f[x]&&(c=f[x][0]);"C"!=c&&(q[x]=c,x&&(t=q[x-1]));f[x]=p(f[x],n,t);"A"!=q[x]&&"C"==c&&(q[x]="C");s(f,x);l&&(l[x]&&(c=l[x][0]),"C"!=c&&(u[x]=c,x&&(t=u[x-1])),l[x]=p(l[x],k,t),"A"!=u[x]&&"C"==c&&(u[x]="C"),s(l,x));r(f,l,n,k,x);r(l,f,k,n,x);var w=f[x],z=l&&l[x],y=w.length,U=l&&z.length;n.x=w[y-2];n.y=w[y-1];n.bx=$(w[y-4])||n.x;n.by=$(w[y-3])||n.y;k.bx=l&&($(z[U-4])||k.x);k.by=l&&($(z[U-3])||k.y);k.x=l&&z[U-2];k.y=l&&z[U-1]}l||(e.curve=d(f));return l?[f,l]:f}function P(a,
b){for(var d=[],e=0,h=a.length;h-2*!b>e;e+=2){var f=[{x:+a[e-2],y:+a[e-1]},{x:+a[e],y:+a[e+1]},{x:+a[e+2],y:+a[e+3]},{x:+a[e+4],y:+a[e+5]}];b?e?h-4==e?f[3]={x:+a[0],y:+a[1]}:h-2==e&&(f[2]={x:+a[0],y:+a[1]},f[3]={x:+a[2],y:+a[3]}):f[0]={x:+a[h-2],y:+a[h-1]}:h-4==e?f[3]=f[2]:e||(f[0]={x:+a[e],y:+a[e+1]});d.push(["C",(-f[0].x+6*f[1].x+f[2].x)/6,(-f[0].y+6*f[1].y+f[2].y)/6,(f[1].x+6*f[2].x-f[3].x)/6,(f[1].y+6*f[2].y-f[3].y)/6,f[2].x,f[2].y])}return d}y=k.prototype;var Q=a.is,C=a._.clone,L="hasOwnProperty",
N=/,?([a-z]),?/gi,$=parseFloat,F=Math,S=F.PI,X=F.min,W=F.max,ma=F.pow,Z=F.abs;M=n(1);var na=n(),ba=n(0,1),V=a._unit2px;a.path=A;a.path.getTotalLength=M;a.path.getPointAtLength=na;a.path.getSubpath=function(a,b,d){if(1E-6>this.getTotalLength(a)-d)return ba(a,b).end;a=ba(a,d,1);return b?ba(a,b).end:a};y.getTotalLength=function(){if(this.node.getTotalLength)return this.node.getTotalLength()};y.getPointAtLength=function(a){return na(this.attr("d"),a)};y.getSubpath=function(b,d){return a.path.getSubpath(this.attr("d"),
b,d)};a._.box=w;a.path.findDotsAtSegment=u;a.path.bezierBBox=p;a.path.isPointInsideBBox=b;a.path.isBBoxIntersect=q;a.path.intersection=function(a,b){return l(a,b)};a.path.intersectionNumber=function(a,b){return l(a,b,1)};a.path.isPointInside=function(a,d,e){var h=r(a);return b(h,d,e)&&1==l(a,[["M",d,e],["H",h.x2+10] ],1)%2};a.path.getBBox=r;a.path.get={path:function(a){return a.attr("path")},circle:function(a){a=V(a);return x(a.cx,a.cy,a.r)},ellipse:function(a){a=V(a);return x(a.cx||0,a.cy||0,a.rx,
a.ry)},rect:function(a){a=V(a);return s(a.x||0,a.y||0,a.width,a.height,a.rx,a.ry)},image:function(a){a=V(a);return s(a.x||0,a.y||0,a.width,a.height)},line:function(a){return"M"+[a.attr("x1")||0,a.attr("y1")||0,a.attr("x2"),a.attr("y2")]},polyline:function(a){return"M"+a.attr("points")},polygon:function(a){return"M"+a.attr("points")+"z"},deflt:function(a){a=a.node.getBBox();return s(a.x,a.y,a.width,a.height)}};a.path.toRelative=function(b){var e=A(b),h=String.prototype.toLowerCase;if(e.rel)return d(e.rel);
a.is(b,"array")&&a.is(b&&b[0],"array")||(b=a.parsePathString(b));var f=[],l=0,n=0,k=0,p=0,s=0;"M"==b[0][0]&&(l=b[0][1],n=b[0][2],k=l,p=n,s++,f.push(["M",l,n]));for(var r=b.length;s<r;s++){var q=f[s]=[],x=b[s];if(x[0]!=h.call(x[0]))switch(q[0]=h.call(x[0]),q[0]){case "a":q[1]=x[1];q[2]=x[2];q[3]=x[3];q[4]=x[4];q[5]=x[5];q[6]=+(x[6]-l).toFixed(3);q[7]=+(x[7]-n).toFixed(3);break;case "v":q[1]=+(x[1]-n).toFixed(3);break;case "m":k=x[1],p=x[2];default:for(var c=1,t=x.length;c<t;c++)q[c]=+(x[c]-(c%2?l:
n)).toFixed(3)}else for(f[s]=[],"m"==x[0]&&(k=x[1]+l,p=x[2]+n),q=0,c=x.length;q<c;q++)f[s][q]=x[q];x=f[s].length;switch(f[s][0]){case "z":l=k;n=p;break;case "h":l+=+f[s][x-1];break;case "v":n+=+f[s][x-1];break;default:l+=+f[s][x-2],n+=+f[s][x-1]}}f.toString=z;e.rel=d(f);return f};a.path.toAbsolute=G;a.path.toCubic=I;a.path.map=function(a,b){if(!b)return a;var d,e,h,f,l,n,k;a=I(a);h=0;for(l=a.length;h<l;h++)for(k=a[h],f=1,n=k.length;f<n;f+=2)d=b.x(k[f],k[f+1]),e=b.y(k[f],k[f+1]),k[f]=d,k[f+1]=e;return a};
a.path.toString=z;a.path.clone=d});C.plugin(function(a,v,y,C){var A=Math.max,w=Math.min,z=function(a){this.items=[];this.bindings={};this.length=0;this.type="set";if(a)for(var f=0,n=a.length;f<n;f++)a[f]&&(this[this.items.length]=this.items[this.items.length]=a[f],this.length++)};v=z.prototype;v.push=function(){for(var a,f,n=0,k=arguments.length;n<k;n++)if(a=arguments[n])f=this.items.length,this[f]=this.items[f]=a,this.length++;return this};v.pop=function(){this.length&&delete this[this.length--];
return this.items.pop()};v.forEach=function(a,f){for(var n=0,k=this.items.length;n<k&&!1!==a.call(f,this.items[n],n);n++);return this};v.animate=function(d,f,n,u){"function"!=typeof n||n.length||(u=n,n=L.linear);d instanceof a._.Animation&&(u=d.callback,n=d.easing,f=n.dur,d=d.attr);var p=arguments;if(a.is(d,"array")&&a.is(p[p.length-1],"array"))var b=!0;var q,e=function(){q?this.b=q:q=this.b},l=0,r=u&&function(){l++==this.length&&u.call(this)};return this.forEach(function(a,l){k.once("snap.animcreated."+
a.id,e);b?p[l]&&a.animate.apply(a,p[l]):a.animate(d,f,n,r)})};v.remove=function(){for(;this.length;)this.pop().remove();return this};v.bind=function(a,f,k){var u={};if("function"==typeof f)this.bindings[a]=f;else{var p=k||a;this.bindings[a]=function(a){u[p]=a;f.attr(u)}}return this};v.attr=function(a){var f={},k;for(k in a)if(this.bindings[k])this.bindings[k](a[k]);else f[k]=a[k];a=0;for(k=this.items.length;a<k;a++)this.items[a].attr(f);return this};v.clear=function(){for(;this.length;)this.pop()};
v.splice=function(a,f,k){a=0>a?A(this.length+a,0):a;f=A(0,w(this.length-a,f));var u=[],p=[],b=[],q;for(q=2;q<arguments.length;q++)b.push(arguments[q]);for(q=0;q<f;q++)p.push(this[a+q]);for(;q<this.length-a;q++)u.push(this[a+q]);var e=b.length;for(q=0;q<e+u.length;q++)this.items[a+q]=this[a+q]=q<e?b[q]:u[q-e];for(q=this.items.length=this.length-=f-e;this[q];)delete this[q++];return new z(p)};v.exclude=function(a){for(var f=0,k=this.length;f<k;f++)if(this[f]==a)return this.splice(f,1),!0;return!1};
v.insertAfter=function(a){for(var f=this.items.length;f--;)this.items[f].insertAfter(a);return this};v.getBBox=function(){for(var a=[],f=[],k=[],u=[],p=this.items.length;p--;)if(!this.items[p].removed){var b=this.items[p].getBBox();a.push(b.x);f.push(b.y);k.push(b.x+b.width);u.push(b.y+b.height)}a=w.apply(0,a);f=w.apply(0,f);k=A.apply(0,k);u=A.apply(0,u);return{x:a,y:f,x2:k,y2:u,width:k-a,height:u-f,cx:a+(k-a)/2,cy:f+(u-f)/2}};v.clone=function(a){a=new z;for(var f=0,k=this.items.length;f<k;f++)a.push(this.items[f].clone());
return a};v.toString=function(){return"Snap\u2018s set"};v.type="set";a.set=function(){var a=new z;arguments.length&&a.push.apply(a,Array.prototype.slice.call(arguments,0));return a}});C.plugin(function(a,v,y,C){function A(a){var b=a[0];switch(b.toLowerCase()){case "t":return[b,0,0];case "m":return[b,1,0,0,1,0,0];case "r":return 4==a.length?[b,0,a[2],a[3] ]:[b,0];case "s":return 5==a.length?[b,1,1,a[3],a[4] ]:3==a.length?[b,1,1]:[b,1]}}function w(b,d,f){d=q(d).replace(/\.{3}|\u2026/g,b);b=a.parseTransformString(b)||
[];d=a.parseTransformString(d)||[];for(var k=Math.max(b.length,d.length),p=[],v=[],h=0,w,z,y,I;h<k;h++){y=b[h]||A(d[h]);I=d[h]||A(y);if(y[0]!=I[0]||"r"==y[0].toLowerCase()&&(y[2]!=I[2]||y[3]!=I[3])||"s"==y[0].toLowerCase()&&(y[3]!=I[3]||y[4]!=I[4])){b=a._.transform2matrix(b,f());d=a._.transform2matrix(d,f());p=[["m",b.a,b.b,b.c,b.d,b.e,b.f] ];v=[["m",d.a,d.b,d.c,d.d,d.e,d.f] ];break}p[h]=[];v[h]=[];w=0;for(z=Math.max(y.length,I.length);w<z;w++)w in y&&(p[h][w]=y[w]),w in I&&(v[h][w]=I[w])}return{from:u(p),
to:u(v),f:n(p)}}function z(a){return a}function d(a){return function(b){return+b.toFixed(3)+a}}function f(b){return a.rgb(b[0],b[1],b[2])}function n(a){var b=0,d,f,k,n,h,p,q=[];d=0;for(f=a.length;d<f;d++){h="[";p=['"'+a[d][0]+'"'];k=1;for(n=a[d].length;k<n;k++)p[k]="val["+b++ +"]";h+=p+"]";q[d]=h}return Function("val","return Snap.path.toString.call(["+q+"])")}function u(a){for(var b=[],d=0,f=a.length;d<f;d++)for(var k=1,n=a[d].length;k<n;k++)b.push(a[d][k]);return b}var p={},b=/[a-z]+$/i,q=String;
p.stroke=p.fill="colour";v.prototype.equal=function(a,b){return k("snap.util.equal",this,a,b).firstDefined()};k.on("snap.util.equal",function(e,k){var r,s;r=q(this.attr(e)||"");var x=this;if(r==+r&&k==+k)return{from:+r,to:+k,f:z};if("colour"==p[e])return r=a.color(r),s=a.color(k),{from:[r.r,r.g,r.b,r.opacity],to:[s.r,s.g,s.b,s.opacity],f:f};if("transform"==e||"gradientTransform"==e||"patternTransform"==e)return k instanceof a.Matrix&&(k=k.toTransformString()),a._.rgTransform.test(k)||(k=a._.svgTransform2string(k)),
w(r,k,function(){return x.getBBox(1)});if("d"==e||"path"==e)return r=a.path.toCubic(r,k),{from:u(r[0]),to:u(r[1]),f:n(r[0])};if("points"==e)return r=q(r).split(a._.separator),s=q(k).split(a._.separator),{from:r,to:s,f:function(a){return a}};aUnit=r.match(b);s=q(k).match(b);return aUnit&&aUnit==s?{from:parseFloat(r),to:parseFloat(k),f:d(aUnit)}:{from:this.asPX(e),to:this.asPX(e,k),f:z}})});C.plugin(function(a,v,y,C){var A=v.prototype,w="createTouch"in C.doc;v="click dblclick mousedown mousemove mouseout mouseover mouseup touchstart touchmove touchend touchcancel".split(" ");
var z={mousedown:"touchstart",mousemove:"touchmove",mouseup:"touchend"},d=function(a,b){var d="y"==a?"scrollTop":"scrollLeft",e=b&&b.node?b.node.ownerDocument:C.doc;return e[d in e.documentElement?"documentElement":"body"][d]},f=function(){this.returnValue=!1},n=function(){return this.originalEvent.preventDefault()},u=function(){this.cancelBubble=!0},p=function(){return this.originalEvent.stopPropagation()},b=function(){if(C.doc.addEventListener)return function(a,b,e,f){var k=w&&z[b]?z[b]:b,l=function(k){var l=
d("y",f),q=d("x",f);if(w&&z.hasOwnProperty(b))for(var r=0,u=k.targetTouches&&k.targetTouches.length;r<u;r++)if(k.targetTouches[r].target==a||a.contains(k.targetTouches[r].target)){u=k;k=k.targetTouches[r];k.originalEvent=u;k.preventDefault=n;k.stopPropagation=p;break}return e.call(f,k,k.clientX+q,k.clientY+l)};b!==k&&a.addEventListener(b,l,!1);a.addEventListener(k,l,!1);return function(){b!==k&&a.removeEventListener(b,l,!1);a.removeEventListener(k,l,!1);return!0}};if(C.doc.attachEvent)return function(a,
b,e,h){var k=function(a){a=a||h.node.ownerDocument.window.event;var b=d("y",h),k=d("x",h),k=a.clientX+k,b=a.clientY+b;a.preventDefault=a.preventDefault||f;a.stopPropagation=a.stopPropagation||u;return e.call(h,a,k,b)};a.attachEvent("on"+b,k);return function(){a.detachEvent("on"+b,k);return!0}}}(),q=[],e=function(a){for(var b=a.clientX,e=a.clientY,f=d("y"),l=d("x"),n,p=q.length;p--;){n=q[p];if(w)for(var r=a.touches&&a.touches.length,u;r--;){if(u=a.touches[r],u.identifier==n.el._drag.id||n.el.node.contains(u.target)){b=
u.clientX;e=u.clientY;(a.originalEvent?a.originalEvent:a).preventDefault();break}}else a.preventDefault();b+=l;e+=f;k("snap.drag.move."+n.el.id,n.move_scope||n.el,b-n.el._drag.x,e-n.el._drag.y,b,e,a)}},l=function(b){a.unmousemove(e).unmouseup(l);for(var d=q.length,f;d--;)f=q[d],f.el._drag={},k("snap.drag.end."+f.el.id,f.end_scope||f.start_scope||f.move_scope||f.el,b);q=[]};for(y=v.length;y--;)(function(d){a[d]=A[d]=function(e,f){a.is(e,"function")&&(this.events=this.events||[],this.events.push({name:d,
f:e,unbind:b(this.node||document,d,e,f||this)}));return this};a["un"+d]=A["un"+d]=function(a){for(var b=this.events||[],e=b.length;e--;)if(b[e].name==d&&(b[e].f==a||!a)){b[e].unbind();b.splice(e,1);!b.length&&delete this.events;break}return this}})(v[y]);A.hover=function(a,b,d,e){return this.mouseover(a,d).mouseout(b,e||d)};A.unhover=function(a,b){return this.unmouseover(a).unmouseout(b)};var r=[];A.drag=function(b,d,f,h,n,p){function u(r,v,w){(r.originalEvent||r).preventDefault();this._drag.x=v;
this._drag.y=w;this._drag.id=r.identifier;!q.length&&a.mousemove(e).mouseup(l);q.push({el:this,move_scope:h,start_scope:n,end_scope:p});d&&k.on("snap.drag.start."+this.id,d);b&&k.on("snap.drag.move."+this.id,b);f&&k.on("snap.drag.end."+this.id,f);k("snap.drag.start."+this.id,n||h||this,v,w,r)}if(!arguments.length){var v;return this.drag(function(a,b){this.attr({transform:v+(v?"T":"t")+[a,b]})},function(){v=this.transform().local})}this._drag={};r.push({el:this,start:u});this.mousedown(u);return this};
A.undrag=function(){for(var b=r.length;b--;)r[b].el==this&&(this.unmousedown(r[b].start),r.splice(b,1),k.unbind("snap.drag.*."+this.id));!r.length&&a.unmousemove(e).unmouseup(l);return this}});C.plugin(function(a,v,y,C){y=y.prototype;var A=/^\s*url\((.+)\)/,w=String,z=a._.$;a.filter={};y.filter=function(d){var f=this;"svg"!=f.type&&(f=f.paper);d=a.parse(w(d));var k=a._.id(),u=z("filter");z(u,{id:k,filterUnits:"userSpaceOnUse"});u.appendChild(d.node);f.defs.appendChild(u);return new v(u)};k.on("snap.util.getattr.filter",
function(){k.stop();var d=z(this.node,"filter");if(d)return(d=w(d).match(A))&&a.select(d[1])});k.on("snap.util.attr.filter",function(d){if(d instanceof v&&"filter"==d.type){k.stop();var f=d.node.id;f||(z(d.node,{id:d.id}),f=d.id);z(this.node,{filter:a.url(f)})}d&&"none"!=d||(k.stop(),this.node.removeAttribute("filter"))});a.filter.blur=function(d,f){null==d&&(d=2);return a.format('<feGaussianBlur stdDeviation="{def}"/>',{def:null==f?d:[d,f]})};a.filter.blur.toString=function(){return this()};a.filter.shadow=
function(d,f,k,u,p){"string"==typeof k&&(p=u=k,k=4);"string"!=typeof u&&(p=u,u="#000");null==k&&(k=4);null==p&&(p=1);null==d&&(d=0,f=2);null==f&&(f=d);u=a.color(u||"#000");return a.format('<feGaussianBlur in="SourceAlpha" stdDeviation="{blur}"/><feOffset dx="{dx}" dy="{dy}" result="offsetblur"/><feFlood flood-color="{color}"/><feComposite in2="offsetblur" operator="in"/><feComponentTransfer><feFuncA type="linear" slope="{opacity}"/></feComponentTransfer><feMerge><feMergeNode/><feMergeNode in="SourceGraphic"/></feMerge>',
{color:u,dx:d,dy:f,blur:k,opacity:p})};a.filter.shadow.toString=function(){return this()};a.filter.grayscale=function(d){null==d&&(d=1);return a.format('<feColorMatrix type="matrix" values="{a} {b} {c} 0 0 {d} {e} {f} 0 0 {g} {b} {h} 0 0 0 0 0 1 0"/>',{a:0.2126+0.7874*(1-d),b:0.7152-0.7152*(1-d),c:0.0722-0.0722*(1-d),d:0.2126-0.2126*(1-d),e:0.7152+0.2848*(1-d),f:0.0722-0.0722*(1-d),g:0.2126-0.2126*(1-d),h:0.0722+0.9278*(1-d)})};a.filter.grayscale.toString=function(){return this()};a.filter.sepia=
function(d){null==d&&(d=1);return a.format('<feColorMatrix type="matrix" values="{a} {b} {c} 0 0 {d} {e} {f} 0 0 {g} {h} {i} 0 0 0 0 0 1 0"/>',{a:0.393+0.607*(1-d),b:0.769-0.769*(1-d),c:0.189-0.189*(1-d),d:0.349-0.349*(1-d),e:0.686+0.314*(1-d),f:0.168-0.168*(1-d),g:0.272-0.272*(1-d),h:0.534-0.534*(1-d),i:0.131+0.869*(1-d)})};a.filter.sepia.toString=function(){return this()};a.filter.saturate=function(d){null==d&&(d=1);return a.format('<feColorMatrix type="saturate" values="{amount}"/>',{amount:1-
d})};a.filter.saturate.toString=function(){return this()};a.filter.hueRotate=function(d){return a.format('<feColorMatrix type="hueRotate" values="{angle}"/>',{angle:d||0})};a.filter.hueRotate.toString=function(){return this()};a.filter.invert=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="table" tableValues="{amount} {amount2}"/><feFuncG type="table" tableValues="{amount} {amount2}"/><feFuncB type="table" tableValues="{amount} {amount2}"/></feComponentTransfer>',{amount:d,
amount2:1-d})};a.filter.invert.toString=function(){return this()};a.filter.brightness=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="linear" slope="{amount}"/><feFuncG type="linear" slope="{amount}"/><feFuncB type="linear" slope="{amount}"/></feComponentTransfer>',{amount:d})};a.filter.brightness.toString=function(){return this()};a.filter.contrast=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="linear" slope="{amount}" intercept="{amount2}"/><feFuncG type="linear" slope="{amount}" intercept="{amount2}"/><feFuncB type="linear" slope="{amount}" intercept="{amount2}"/></feComponentTransfer>',
{amount:d,amount2:0.5-d/2})};a.filter.contrast.toString=function(){return this()}});return C});

]]> </script>
</svg>



<div class="markdown"><p>Some methods are implemented by searching this tree in depth first order</p>
</div>

<pre class='language-julia'><code class='language-julia'>gatecount(qft3)</code></pre>
<pre id='var-hash311647' class='code-output documenter-example-output'>Dict{Type{<:AbstractBlock}, Int64} with 2 entries:
  HGate                                  => 3
  ControlBlock{ShiftGate{Float64}, 1, 1} => 3</pre>

<pre class='language-julia'><code class='language-julia'>parameters(qft3)</code></pre>
<pre id='var-hash163956' class='code-output documenter-example-output'>3-element Vector{Float64}:
 1.5707963267948966
 0.7853981633974483
 1.5707963267948966</pre>

<pre class='language-julia'><code class='language-julia'>dispatch(qft3, [1.0, 2.0, 3.0])</code></pre>
<pre id='var-hash115737' class='code-output documenter-example-output'>nqubits: 3
chain
├─ chain
│  ├─ put on (1)
│  │  └─ H
│  ├─ control(2)
│  │  └─ (1,) shift(1.0)
│  └─ control(3)
│     └─ (1,) shift(2.0)
├─ chain
│  ├─ put on (2)
│  │  └─ H
│  └─ control(3)
│     └─ (2,) shift(3.0)
└─ chain
   └─ put on (3)
      └─ H
</pre>


<div class="markdown"><h3>Simplifying the block tree</h3>
</div>


<div class="markdown"><p>If you feel hard to handle two many types, you can compile away non-essential composite blocks like <code>repeat</code> and <code>kron</code> to basic types.</p>
</div>

<pre class='language-julia'><code class='language-julia'>Optimise.to_basictypes(kron(6, 1=&gt;X, (2,3)=&gt;SWAP))</code></pre>
<pre id='var-hash693750' class='code-output documenter-example-output'>nqubits: 6
chain
├─ put on (1)
│  └─ X
└─ put on (2, 3)
   └─ SWAP
</pre>

<pre class='language-julia'><code class='language-julia'>kron(6, 1=&gt;X, (2,3)=&gt;SWAP)</code></pre>
<pre id='var-hash162061' class='code-output documenter-example-output'>nqubits: 6
kron
├─ 1=>X
└─ 2:3=>SWAP</pre>

<pre class='language-julia'><code class='language-julia'>Optimise.to_basictypes(repeat(10, X, (5,2)))</code></pre>
<pre id='var-hash229770' class='code-output documenter-example-output'>nqubits: 10
chain
├─ put on (2)
│  └─ X
└─ put on (5)
   └─ X
</pre>


<div class="markdown"><p>However, these function does not propagate.</p>
</div>

<pre class='language-julia'><code class='language-julia'>Optimise.to_basictypes(chain(repeat(10, X, (5,2))))</code></pre>
<pre id='var-hash427405' class='code-output documenter-example-output'>nqubits: 10
chain
└─ repeat on (2, 5)
   └─ X
</pre>


<div class="markdown"><p>One apply these functions as the simplification rules repeatedly through the tree</p>
</div>

<pre class='language-julia'><code class='language-julia'>Optimise.simplify(chain(repeat(10, X, (5,2))), rules=[Optimise.to_basictypes])</code></pre>
<pre id='var-hash146896' class='code-output documenter-example-output'>nqubits: 10
chain
└─ chain
   ├─ put on (2)
   │  └─ X
   └─ put on (5)
      └─ X
</pre>


<div class="markdown"><p>Apply this rule to the Heisenberg Hamiltonian</p>
</div>

<pre class='language-julia'><code class='language-julia'>Optimise.simplify(EasyBuild.heisenberg(3), rules=[Optimise.to_basictypes])</code></pre>
<pre id='var-hash502633' class='code-output documenter-example-output'>nqubits: 3
+
├─ +
│  ├─ chain
│  │  ├─ put on (1)
│  │  │  └─ X
│  │  └─ put on (2)
│  │     └─ X
│  ├─ chain
│  │  ├─ put on (1)
│  │  │  └─ Y
│  │  └─ put on (2)
│  │     └─ Y
│  └─ chain
│     ├─ put on (1)
│     │  └─ Z
│     └─ put on (2)
│        └─ Z
├─ +
│  ├─ chain
│  │  ├─ put on (2)
│  │  │  └─ X
│  │  └─ put on (3)
│  │     └─ X
│  ├─ chain
│  │  ├─ put on (2)
│  │  │  └─ Y
│  │  └─ put on (3)
│  │     └─ Y
│  └─ chain
│     ├─ put on (2)
│     │  └─ Z
│     └─ put on (3)
│        └─ Z
└─ +
   ├─ chain
   │  ├─ put on (1)
   │  │  └─ X
   │  └─ put on (3)
   │     └─ X
   ├─ chain
   │  ├─ put on (1)
   │  │  └─ Y
   │  └─ put on (3)
   │     └─ Y
   └─ chain
      ├─ put on (1)
      │  └─ Z
      └─ put on (3)
         └─ Z
</pre>


<div class="markdown"><p>You can also remove nested <code>ChainBlock</code> or <code>Add</code> block</p>
</div>

<pre class='language-julia'><code class='language-julia'>Optimise.simplify(EasyBuild.heisenberg(3), rules=[Optimise.to_basictypes, Optimise.eliminate_nested])</code></pre>
<pre id='var-hash155025' class='code-output documenter-example-output'>nqubits: 3
+
├─ chain
│  ├─ put on (1)
│  │  └─ X
│  └─ put on (2)
│     └─ X
├─ chain
│  ├─ put on (1)
│  │  └─ Y
│  └─ put on (2)
│     └─ Y
├─ chain
│  ├─ put on (1)
│  │  └─ Z
│  └─ put on (2)
│     └─ Z
├─ chain
│  ├─ put on (2)
│  │  └─ X
│  └─ put on (3)
│     └─ X
├─ chain
│  ├─ put on (2)
│  │  └─ Y
│  └─ put on (3)
│     └─ Y
├─ chain
│  ├─ put on (2)
│  │  └─ Z
│  └─ put on (3)
│     └─ Z
├─ chain
│  ├─ put on (1)
│  │  └─ X
│  └─ put on (3)
│     └─ X
├─ chain
│  ├─ put on (1)
│  │  └─ Y
│  └─ put on (3)
│     └─ Y
└─ chain
   ├─ put on (1)
   │  └─ Z
   └─ put on (3)
      └─ Z
</pre>


<div class="markdown"><h2>Play with Hamiltonians</h2>
</div>

<pre class='language-julia'><code class='language-julia'>h30 = EasyBuild.heisenberg(30)</code></pre>
<pre id='var-h30' class='code-output documenter-example-output'>nqubits: 30
+
├─ +
│  ├─ repeat on (1, 2)
│  │  └─ X
│  ├─ repeat on (1, 2)
│  │  └─ Y
│  └─ repeat on (1, 2)
│     └─ Z
├─ +
│  ├─ repeat on (2, 3)
│  │  └─ X
│  ├─ repeat on (2, 3)
│  │  └─ Y
│  └─ repeat on (2, 3)
│     └─ Z
├─ +
│  ├─ repeat on (3, 4)
│  │  └─ X
│  ├─ repeat on (3, 4)
│  │  └─ Y
│  └─ repeat on (3, 4)
│     └─ Z
├─ +
│  ├─ repeat on (4, 5)
│  │  └─ X
│  ├─ repeat on (4, 5)
│  │  └─ Y
│  └─ repeat on (4, 5)
│     └─ Z
├─ +
│  ├─ repeat on (5, 6)
│  │  └─ X
│  ├─ repeat on (5, 6)
│  │  └─ Y
│  └─ repeat on (5, 6)
│     └─ Z
├─ +
│  ├─ repeat on (6, 7)
│  │  └─ X
│  ├─ repeat on (6, 7)
│  │  └─ Y
│  └─ repeat on (6, 7)
│     └─ Z
├─ +
│  ├─ repeat on (7, 8)
│  │  └─ X
│  ├─ repeat on (7, 8)
│  │  └─ Y
│  └─ repeat on (7, 8)
│     └─ Z
├─ +
│  ├─ repeat on (8, 9)
│  │  └─ X
│  ├─ repeat on (8, 9)
│  │  └─ Y
│  └─ repeat on (8, 9)
│     └─ Z
├─ +
│  ├─ repeat on (9, 10)
│  │  └─ X
│  ├─ repeat on (9, 10)
│  │  └─ Y
│  └─ repeat on (9, 10)
│     └─ Z
├─ +
│  ├─ repeat on (10, 11)
│  │  └─ X
│  ├─ repeat on (10, 11)
│  │  └─ Y
│  └─ repeat on (10, 11)
│     └─ Z
├─ +
│  ├─ repeat on (11, 12)
│  │  └─ X
│  ├─ repeat on (11, 12)
│  │  └─ Y
│  └─ repeat on (11, 12)
│     └─ Z
├─ +
│  ├─ repeat on (12, 13)
│  │  └─ X
│  ├─ repeat on (12, 13)
│  │  └─ Y
│  └─ repeat on (12, 13)
│     └─ Z
├─ +
│  ├─ repeat on (13, 14)
│  │  └─ X
│  ├─ repeat on (13, 14)
│  │  └─ Y
│  └─ repeat on (13, 14)
│     └─ Z
├─ +
│  ├─ repeat on (14, 15)
│  │  └─ X
│  ├─ repeat on (14, 15)
│  │  └─ Y
│  └─ repeat on (14, 15)
│     └─ Z
├─ +
│  ├─ repeat on (15, 16)
│  │  └─ X
│  ├─ repeat on (15, 16)
│  │  └─ Y
│  └─ repeat on (15, 16)
│     └─ Z
├─ +
│  ├─ repeat on (16, 17)
│  │  └─ X
│  ├─ repeat on (16, 17)
│  │  └─ Y
│  └─ repeat on (16, 17)
│     └─ Z
├─ +
│  ├─ repeat on (17, 18)
│  │  └─ X
│  ├─ repeat on (17, 18)
│  │  └─ Y
│  └─ repeat on (17, 18)
│     └─ Z
├─ +
│  ├─ repeat on (18, 19)
│  │  └─ X
│  ├─ repeat on (18, 19)
│  │  └─ Y
│  └─ repeat on (18, 19)
│     └─ Z
├─ +
│  ├─ repeat on (19, 20)
│  │  └─ X
│  ├─ repeat on (19, 20)
│  │  └─ Y
│  └─ repeat on (19, 20)
│     └─ Z
├─ +
│  ├─ repeat on (20, 21)
│  │  └─ X
│  ├─ repeat on (20, 21)
│  │  └─ Y
│  └─ repeat on (20, 21)
│     └─ Z
├─ +
│  ├─ repeat on (21, 22)
│  │  └─ X
│  ├─ repeat on (21, 22)
│  │  └─ Y
│  └─ repeat on (21, 22)
│     └─ Z
├─ +
│  ├─ repeat on (22, 23)
│  │  └─ X
│  ├─ repeat on (22, 23)
│  │  └─ Y
│  └─ repeat on (22, 23)
│     └─ Z
├─ +
│  ├─ repeat on (23, 24)
│  │  └─ X
│  ├─ repeat on (23, 24)
│  │  └─ Y
│  └─ repeat on (23, 24)
│     └─ Z
├─ +
│  ├─ repeat on (24, 25)
│  │  └─ X
│  ├─ repeat on (24, 25)
│  │  └─ Y
│  └─ repeat on (24, 25)
│     └─ Z
├─ +
│  ├─ repeat on (25, 26)
│  │  └─ X
│  ├─ repeat on (25, 26)
│  │  └─ Y
│  └─ repeat on (25, 26)
│     └─ Z
├─ +
│  ├─ repeat on (26, 27)
│  │  └─ X
│  ├─ repeat on (26, 27)
│  │  └─ Y
│  └─ repeat on (26, 27)
│     └─ Z
├─ +
│  ├─ repeat on (27, 28)
│  │  └─ X
│  ├─ repeat on (27, 28)
│  │  └─ Y
│  └─ repeat on (27, 28)
│     └─ Z
├─ +
│  ├─ repeat on (28, 29)
│  │  └─ X
│  ├─ repeat on (28, 29)
│  │  └─ Y
│  └─ repeat on (28, 29)
│     └─ Z
├─ +
│  ├─ repeat on (29, 30)
│  │  └─ X
│  ├─ repeat on (29, 30)
│  │  └─ Y
│  └─ repeat on (29, 30)
│     └─ Z
└─ +
   ├─ repeat on (1, 30)
   │  └─ X
   ├─ repeat on (1, 30)
   │  └─ Y
   └─ repeat on (1, 30)
      └─ Z
</pre>


<div class="markdown"><h3>Block indexing</h3>
</div>

<pre class='language-julia'><code class='language-julia'># j = DitStr{2}(rand(0:1, 30))  # programming way
# j = dit"0000000000_0000100000_0000000000;2"  # dit string
j = bit"0000000000_0000100000_0000000000"</code></pre>
<pre id='var-j' class='code-output documenter-example-output'>000000000000001000000000000000 ₍₂₎</pre>

<pre class='language-julia'><code class='language-julia'>h30[j, j]</code></pre>
<pre id='var-hash107338' class='code-output documenter-example-output'>26.0 + 0.0im</pre>

<pre class='language-julia'><code class='language-julia'>h30[:, j]</code></pre>
<pre id='var-hash704840' class='code-output documenter-example-output'>YaoBlocks.EntryTable{BitBasis.DitStr{2, 30, Int64}, ComplexF64}:
  000000000000000100000000000000 ₍₂₎   2.0 + 0.0im
  000000000000001000000000000000 ₍₂₎   26.0 + 0.0im
  000000000000010000000000000000 ₍₂₎   2.0 + 0.0im
</pre>

<pre class='language-julia'><code class='language-julia'>h30[:, h30[:, j]]</code></pre>
<pre id='var-hash161055' class='code-output documenter-example-output'>YaoBlocks.EntryTable{BitBasis.DitStr{2, 30, Int64}, ComplexF64}:
  000000000000000010000000000000 ₍₂₎   4.0 + 0.0im
  000000000000000100000000000000 ₍₂₎   104.0 + 0.0im
  000000000000001000000000000000 ₍₂₎   684.0 + 0.0im
  000000000000010000000000000000 ₍₂₎   104.0 + 0.0im
  000000000000100000000000000000 ₍₂₎   4.0 + 0.0im
</pre>

<pre class='language-julia'><code class='language-julia'>let
    gadget = @bind num_propagate Slider(0:20; show_value=true)
    md"""num_propagate $gadget"""
end</code></pre>
<div class="markdown"><p>num_propagate <bond def="num_propagate" unique_id="V0DSPnbjez2x"><input type='range' min='1' max='21' value='1'><script>
					const input_el = currentScript.previousElementSibling
					const output_el = currentScript.nextElementSibling
					const displays = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"]
					
					input_el.addEventListener("input", () => {
						output_el.value = displays[input_el.valueAsNumber - 1]
					})
					</script><output style='
						font-family: system-ui;
    					font-size: 15px;
    					margin-left: 3px;
    					transform: translateY(-4px);
    					display: inline-block;'>0</output></bond></p>
</div>

<pre class='language-julia'><code class='language-julia'>let J = EntryTable([j], [1.0+0im])
    for i=1:num_propagate
        J = h30[:,J]
    end
    J
end</code></pre>
<pre id='var-hash295876' class='code-output documenter-example-output'>YaoBlocks.EntryTable{BitBasis.DitStr{2, 30, Int64}, ComplexF64}:
  000000000000001000000000000000 ₍₂₎   1.0 + 0.0im
</pre>


<div class="markdown"><h3>The ground state</h3>
</div>

<pre class='language-julia'><code class='language-julia'>@time mat(h15)</code></pre>
<pre id='var-hash152416' class='code-output documenter-example-output'>32768×32768 SparseMatrixCSC{ComplexF64, Int64} with 278528 stored entries:
⠻⣦⣄⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠳⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⢙⡻⣮⡳⡄⠀⡀⠀⠀⠀⠀⠀⠀⠀⠈⠳⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠙⠮⢻⣶⡄⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠈⠳⣄⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠠⣄⠉⠻⣦⣀⠙⠦⠀⠀⠀⣄⠀⠀⠀⠀⠈⠳⣄⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠈⠳⣄⠘⢿⣷⡰⣄⠀⠀⠈⠳⣄⠀⠀⠀⠀⠈⠳⣄⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠈⠃⠐⢮⡻⣮⣇⠀⠀⠀⠈⠳⣄⠀⠀⠀⠀⠈⠳⣄⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠙⠿⣧⠀⠀⠀⠀⠈⠳⣄⠀⠀⠀⠀⠈⠳⣄
⠙⢦⡀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⢻⣶⣄⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠙⢦⡀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⢹⡻⣮⡳⠄⢠⡀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠙⠎⢿⣷⡄⠙⢦⡀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠙⠀⠀⠀⠲⣄⠉⠻⣦⣀⠙⠂⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠈⠳⣄⠘⠿⣧⡲⣄⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⡀⠀⠀⠀⠀⠀⠀⠀⠈⠀⠘⢮⡻⣮⣅⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠙⠻⣦</pre>

<pre class='language-julia'><code class='language-julia'>@time KrylovKit.eigsolve(mat(h15), statevec(rand_state(15)), 1, :SR)</code></pre>
<pre id='var-hash992904' class='code-output documenter-example-output'>([-26.134670289862953], Vector{ComplexF64}[[3.1614595678894303e-19 - 3.746214512625886e-18im, -9.93407964144914e-19 + 6.08848870649073e-18im, 3.57868006165088e-18 + 1.6344235984169927e-19im, -1.9716634488992694e-19 - 2.2753670338007436e-18im, -2.8204591823381953e-18 + 6.114048267666856e-19im, 2.1614916679238586e-18 + 7.037495268733614e-18im, 1.993376267320287e-18 - 9.025936766281825e-19im, -2.3105587014319464e-18 + 4.70409223951071e-19im, -2.623418839018243e-18 - 9.299381101751635e-19im, 2.7735108643057407e-18 - 2.913196654947031e-18im  …  6.915257068269409e-18 - 5.038730326586602e-18im, -5.264920644357111e-19 + 1.4696059132096758e-18im, -2.2597246220453763e-18 - 4.382331397148873e-18im, -3.7929251569545465e-19 + 1.7182095963415168e-18im, 8.000312726245381e-18 - 1.0431588007365672e-18im, 9.347326914902176e-18 + 1.8080149282637537e-19im, -1.5576694376799531e-18 - 4.543662235187159e-18im, -2.3724480417208914e-18 - 5.1842552035075844e-18im, -1.8153975848552832e-19 + 4.30043909766949e-18im, 4.943091443268305e-18 - 8.88905887450899e-18im]], ConvergenceInfo: one converged value after 4 iterations and 66 applications of the linear map;
norms of residuals are given by (3.87215331946589e-15,).
)</pre>


<div class="markdown"><h3>Time Evolution</h3>
</div>

<pre class='language-julia'><code class='language-julia'>h15 = EasyBuild.heisenberg(15);</code></pre>


<pre class='language-julia'><code class='language-julia'>@time let
    reg = rand_state(15)
    apply(reg, time_evolve(h15, 0.5))
end</code></pre>
<pre id='var-hash152053' class='code-output documenter-example-output'>ArrayReg{2, ComplexF64, Array...}
    active qubits: 15/15
    nlevel: 2</pre>


<div class="markdown"><p>For larger Hamiltonian, we prefer storing the Hamiltonian matrix in the memory instead of computing the terms on the fly.</p>
</div>

<pre class='language-julia'><code class='language-julia'>cache(h15)</code></pre>
<pre id='var-hash150597' class='code-output documenter-example-output'>nqubits: 15
[cached] +
   ├─ +
   │  ├─ repeat on (1, 2)
   │  │  └─ X
   │  ├─ repeat on (1, 2)
   │  │  └─ Y
   │  └─ repeat on (1, 2)
   │     └─ Z
   ├─ +
   │  ├─ repeat on (2, 3)
   │  │  └─ X
   │  ├─ repeat on (2, 3)
   │  │  └─ Y
   │  └─ repeat on (2, 3)
   │     └─ Z
   ├─ +
   │  ├─ repeat on (3, 4)
   │  │  └─ X
   │  ├─ repeat on (3, 4)
   │  │  └─ Y
   │  └─ repeat on (3, 4)
   │     └─ Z
   ├─ +
   │  ├─ repeat on (4, 5)
   │  │  └─ X
   │  ├─ repeat on (4, 5)
   │  │  └─ Y
   │  └─ repeat on (4, 5)
   │     └─ Z
   ├─ +
   │  ├─ repeat on (5, 6)
   │  │  └─ X
   │  ├─ repeat on (5, 6)
   │  │  └─ Y
   │  └─ repeat on (5, 6)
   │     └─ Z
   ├─ +
   │  ├─ repeat on (6, 7)
   │  │  └─ X
   │  ├─ repeat on (6, 7)
   │  │  └─ Y
   │  └─ repeat on (6, 7)
   │     └─ Z
   ├─ +
   │  ├─ repeat on (7, 8)
   │  │  └─ X
   │  ├─ repeat on (7, 8)
   │  │  └─ Y
   │  └─ repeat on (7, 8)
   │     └─ Z
   ├─ +
   │  ├─ repeat on (8, 9)
   │  │  └─ X
   │  ├─ repeat on (8, 9)
   │  │  └─ Y
   │  └─ repeat on (8, 9)
   │     └─ Z
   ├─ +
   │  ├─ repeat on (9, 10)
   │  │  └─ X
   │  ├─ repeat on (9, 10)
   │  │  └─ Y
   │  └─ repeat on (9, 10)
   │     └─ Z
   ├─ +
   │  ├─ repeat on (10, 11)
   │  │  └─ X
   │  ├─ repeat on (10, 11)
   │  │  └─ Y
   │  └─ repeat on (10, 11)
   │     └─ Z
   ├─ +
   │  ├─ repeat on (11, 12)
   │  │  └─ X
   │  ├─ repeat on (11, 12)
   │  │  └─ Y
   │  └─ repeat on (11, 12)
   │     └─ Z
   ├─ +
   │  ├─ repeat on (12, 13)
   │  │  └─ X
   │  ├─ repeat on (12, 13)
   │  │  └─ Y
   │  └─ repeat on (12, 13)
   │     └─ Z
   ├─ +
   │  ├─ repeat on (13, 14)
   │  │  └─ X
   │  ├─ repeat on (13, 14)
   │  │  └─ Y
   │  └─ repeat on (13, 14)
   │     └─ Z
   ├─ +
   │  ├─ repeat on (14, 15)
   │  │  └─ X
   │  ├─ repeat on (14, 15)
   │  │  └─ Y
   │  └─ repeat on (14, 15)
   │     └─ Z
   └─ +
      ├─ repeat on (1, 15)
      │  └─ X
      ├─ repeat on (1, 15)
      │  └─ Y
      └─ repeat on (1, 15)
         └─ Z
</pre>

<pre class='language-julia'><code class='language-julia'>@time let
    reg = rand_state(15)
    apply(reg, time_evolve(cache(h15), 0.5))
end</code></pre>
<pre id='var-hash106915' class='code-output documenter-example-output'>ArrayReg{2, ComplexF64, Array...}
    active qubits: 15/15
    nlevel: 2</pre>


<div class="markdown"><h2>Eigen basis and measurement</h2>
</div>

<pre class='language-julia'><code class='language-julia'>YaoBlocks.eigenbasis(kron(4, 4=&gt;X, 1=&gt;X))</code></pre>
<pre id='var-hash174524' class='code-output documenter-example-output'>(nqubits: 4
kron
├─ 1=>Z
└─ 4=>Z, nqubits: 4
kron
├─ 1=>H
└─ 4=>H)</pre>

<pre class='language-julia'><code class='language-julia'>let
    reg = rand_state(4)
    res = measure!(reg)
    print_table(reg)
    res
end</code></pre>
<pre id='var-hash104709' class='code-output documenter-example-output'>0111 ₍₂₎</pre>

<pre class='language-julia'><code class='language-julia'>let
    reg = rand_state(4)
    res = measure!(kron(4, 4=&gt;X, 1=&gt;X), reg)
    print_table(apply(reg, kron(4, 4=&gt;H, 1=&gt;H)))  # print on the eigenbasis
    res  # the state projected to that eigenbasis
end</code></pre>
<pre id='var-hash293010' class='code-output documenter-example-output'>-1.0 + 0.0im</pre>


<div class="markdown"><h2>Extend the block system</h2>
</div>


<div class="markdown"><p><strong>General purposed</strong></p>
<ul>
<li><p><code>YaoAPI.mat&#40;T, block&#41; -&gt; AbstractMatrix</code>, matrix representation</p>
</li>
<li><p><code>YaoBlocks.unsafe_apply&#33;&#40;reg, block&#41;</code>, allocation free method to apply a gate &#40;optional&#41;</p>
</li>
<li><p><code>YaoBlocks.print_block&#40;io, block&#41;</code>, better printing &#40;optional&#41;</p>
</li>
</ul>
<p><strong>Composite blocks</strong></p>
<ul>
<li><p><code>YaoAPI.occupied_locs&#40;block&#41; -&gt; tuple</code></p>
</li>
<li><p><code>YaoAPI.subblocks&#40;block&#41; -&gt; tuple or vector</code></p>
</li>
<li><p><code>YaoAPI.chsubblocks&#40;block, subblocks&#41;</code></p>
</li>
</ul>
<p><strong>Parameter management API</strong></p>
<ul>
<li><p><code>setiparams&#33;&#40;block, parameters&#41;</code>, setter of intrinsic parameters</p>
</li>
<li><p><code>getiparams&#40;block&#41; -&gt; tuple</code>, getter of intrinsic parameters</p>
</li>
</ul>
</div>


<div class="markdown"><p>An exmaple: <a href="https://github.com/QuantumBFS/Yao.jl/blob/master/src/EasyBuild/block_extension/FSimGate.jl">https://github.com/QuantumBFS/Yao.jl/blob/master/src/EasyBuild/block_extension/FSimGate.jl</a></p>
</div>


<div class="markdown"><h2>News</h2>
<ul>
<li><p>Aqua test: <a href="https://github.com/QuantumBFS/Yao.jl/issues/381">https://github.com/QuantumBFS/Yao.jl/issues/381</a></p>
</li>
<li><p>Next time: Bloqade <a href="https://github.com/QuEraComputing/Bloqade.jl">https://github.com/QuEraComputing/Bloqade.jl</a></p>
</li>
<li><p>Announce of Bounty issue: <a href="https://github.com/QuantumBFS/Yao.jl/issues/403">https://github.com/QuantumBFS/Yao.jl/issues/403</a></p>
</li>
</ul>
</div>
<div class='manifest-versions'>
<p>Built with Julia 1.8.0-rc1 and</p>
BenchmarkTools 1.3.1<br>
KrylovKit 0.5.4<br>
PlutoUI 0.7.39<br>
Yao 0.8.2<br>
YaoPlots 0.7.5
</div>

<!-- PlutoStaticHTML.End -->
```
