const TRANSPARENT = 'rgba(0,0,0,0)';

const mapID = 'map-container';
const specID = 'spectrum-container';
const channelID = 'channel';
const chLabelID = 'channelLabel';

let cube = null;

class SpectralCube {
	constructor(wave, flux, metadata) {
		this.wave = wave;
		this.flux = flux;
		this.metadata = metadata;
		[this.nz, this.ny, this.nx] = this.metadata.shape;
		this.channel = 0;
		this.focusX = 0;
		this.focusY = 0;

		this.mapContainer = document.getElementById(mapID);
		this.channelSlider = document.getElementById(channelID);
		this.chLabel = document.getElementById(chLabelID);

		this.update();
		this.initEventHandlers();
	}

	updateFocus(x, y) {
		const newx = Math.max(0, Math.min(this.nx-1,Math.round(x)));
		const newy = Math.max(0, Math.min(this.ny-1,Math.round(y)));
		this.focusX = newx;
		this.focusY = newy;
		this.update();
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

		// console.log(map);
		return map;
	}


	update() {
		this.updateMap();
		this.updateSpectrum();
	}


	updateMap() {
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

		const heatmap = {z:this.getMap(this.channel),type:'heatmap'};
		const focusBox = {
			type: 'scatter', mode: 'markers',
			x: [this.focusX], y: [this.focusY],
			marker: {
				symbol: 'x-thin', color: '#ffffff',
				line: { color: '#ffffff', width: 1.0 },
			},
			hoverinfo: 'skip', showlegend: false
		};

		Plotly.react(mapID, [heatmap, focusBox], layout, {responsive:true});
		this.chLabel.textContent = this.wave[this.channel].toPrecision(5).toString();
	}


	updateSpectrum(){
		let spax = this.getSpaxSpec(this.focusX,this.focusY);
		var plot = {
			type: 'scatter', mode: 'lines',
			x: this.wave,
			y: spax,
		};
		var layout = {
			xaxis: {title: 'Wavelength'}, yaxis: {title: 'Flux'}
		};
		Plotly.react(specID, [plot], layout, {responsive:true});
	}


	onMapClick(event) {
		if (event?.points?.length) {
			this.updateFocus(event.points[0].x, event.points[0].y);
		}
	}


	onChannelSlide(event) {
		let ch = Number(event.target.value);
		this.channel = Math.max(0, Math.min(this.nz-1, ch));
		this.updateMap();
	}


	onKeyDown(event) {
		if (['INPUT', 'TEXTAREA', 'SELECT'].includes(document.activeElement?.tagName)) { return; }

		if (event.key === 'ArrowUp') { this.updateFocus(this.focusX, this.focusY+1); }
		else if (event.key === 'ArrowDown') { this.updateFocus(this.focusX, this.focusY-1); }
		else if (event.key === 'ArrowLeft') { this.updateFocus(this.focusX-1, this.focusY); }
		else if (event.key === 'ArrowRight') { this.updateFocus(this.focusX+1, this.focusY); }
	}


	initEventHandlers() {
		this.mapContainer.on('plotly_click', (e) => {this.onMapClick(e);});
		this.channelSlider.max = this.nz-1;
		this.channelSlider.addEventListener('input', (e) => {this.onChannelSlide(e);});
		window.addEventListener('keydown', (e) => {this.onKeyDown(e);});
	}
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
}


initialize();
