const TRANSPARENT = 'rgba(0,0,0,0)';

const mapID = 'map-container';
const specID = 'spectrum-container';

let cube = null;

class SpectralCube {
	constructor(wave, flux, metadata) {
		this.wave = wave;
		this.flux = flux;
		this.metadata = metadata;
		[this.nz, this.ny, this.nx] = this.metadata.shape;
		this.focusX = 0;
		this.focusY = 0;
	}

	updateFocus(x, y) {
		console.log('updateFocus');
		console.log(x);
		console.log(y);
		this.focusX = x;
		this.focusY = y;
		updateMap(10);
		updateSpectrum();
	}

	getSpaxSpec(x, y) {
		const spec = new Float32Array(this.nz);
		for (let w = 0; w < this.nz; w++) {
			const idx = (w * this.ny + y) * this.nx + x;
			spec[w] = this.flux[idx];
		}
		return spec;
	}

	getMap(w) {
		const s = w * this.nx * this.ny;
		const slice = this.flux.subarray(s, s + this.nx * this.ny);

		let map = new Array(this.ny);
		for (let y = 0; y < this.ny; y++) {
			map[y] = Array.from(slice.subarray(y * this.nx, (y+1) * this.nx));
		}

		map = map.map(row => row.map(val => Math.log10(val)));

		console.log(map);
		return map;
	}
}


function updateMap(w) {
	if (!cube) { return; }

	var layout = {
		title: {text: 'Cube Viewer', subtitle: {text: 'Subtitle'}},
		annotations: [],
		paper_bgcolor: TRANSPARENT,
		plot_bgcolor: TRANSPARENT,
		font: {color: 'black', size:16},
		xaxis: {title: {text:'X (px)'}, gridcolor: 'black', autoresize: true},
		yaxis: {title: {text:'Y (px)'}, gridcolor: 'black', autoresize: true},
		clickmode: 'event+select'
	};

	const heatmap = {z:cube.getMap(w),type:'heatmap'};
	const focusBox = {
		type: 'scatter', mode: 'markers',
		x: [cube.focusX], y: [cube.focusY],
		marker: {
			symbol: 'x-thin', color: '#ffffff',
			line: { color: '#ffffff', width: 1.0 },
		},
		hoverinfo: 'skip', showlegend: false
	};

	Plotly.react(mapID, [heatmap, focusBox], layout, {responsive:true});
}


function updateSpectrum(){
	if (!cube) { return; }

	console.log('updateSpectrum');
	let spax = cube.getSpaxSpec(cube.focusX,cube.focusY);
	var plot = {
		type: 'scatter', mode: 'lines',
		x: cube.wave,
		y: spax,
	};
	var layout = {
		xaxis: {title: 'Wavelength'}, yaxis: {title: 'Flux'}
	};
	Plotly.react(specID, [plot], layout, {responsive:true});
}


function handleClick(e) {
	if (!cube) { return; }

	if (e?.points?.length) {
		x = Math.max(0, Math.min(cube.nx-1,Math.round(e.points[0].x)));
		y = Math.max(0, Math.min(cube.ny-1,Math.round(e.points[0].y)));
		cube.updateFocus(x, y);
	}
}


function initializeMap() {
	updateMap(10);
}


function initializeSpectrum() {
	updateSpectrum()
}


async function loadSpectralCube() {
    const metadata = await fetch('metadata.json')
        .then(response => response.json());

    const [wavelengthBuffer, fluxBuffer] = await Promise.all([
        fetch('./bin/wavelength.bin')
            .then(response => response.arrayBuffer()),

        fetch('./bin/flux.bin')
            .then(response => response.arrayBuffer())
    ]);

    const wave = new Float32Array(wavelengthBuffer);
    const flux = new Float32Array(fluxBuffer);

    cube = new SpectralCube(wave, flux, metadata)
}


async function initialize() {
	await loadSpectralCube();
	console.log(cube);
	initializeMap();
	initializeSpectrum();

	const mapContainer = document.getElementById(mapID);
	mapContainer.on('plotly_click', handleClick);
}


initialize();
