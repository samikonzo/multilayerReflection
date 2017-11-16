'use strict'
var l = console.log;

var base_lambda = 1,
	layer0 = { n : 1 },
	layer1 = { n : 2.2 ,  h: 0.45454545454545454545454545454545}, 
	layer2 = { n : 1.45 , h: 0.68965517241379310344827586206897}, 
	layer3 = { n : 2.2 ,  h: 0.45454545454545454545454545454545}, 
	layer4 = { n : 1.52},

	layer_n2 = { n : 1.45, h: 0.15/1.45};

layer_n2 = { n : 2.2, h: 0.15/2.2};


var options = {
	layers: [
		layer0,
		layer1,
		layer2,
		layer3,
		layer4,
	],

	lambdaStart : 1,
	lambdaEnd : 1,
	lambdaDelta : .1,
}


options = {
	layers: [
		layer0,
		layer_n2,
		layer4,
	],

	lambdaStart : .4,
	lambdaEnd : .75,
	lambdaDelta : .05,
}


calculateReflection(options)

function calculateReflection(options){
	//l(options)

	var layers = options.layers,
		layerCount = layers.length - 1, // count
		m = layerCount, // temp symbol
		lambdaStart = options.lambdaStart,
		lambdaEnd = options.lambdaEnd,
		lambdaDelta = options.lambdaDelta,

		reflection = {}, // r
		deltaPhase = {}, // ∆
		deg_alpha = {}, //α  
		deg_beta = {}, // β
		deg_gama = {}, // γ
		Ri = {};

	const PI = Math.PI;;	

	//calcualte reflection for all wavelength
	for(var lambda = lambdaStart; lambda <= lambdaEnd; lambda += lambdaDelta){
		l(lambda + 'мкм')

		//zeroing out global variables
		reflection = {}, // r
		deltaPhase = {}, // ∆
		deg_alpha = {}, //α  
		deg_beta = {}, // β
		deg_gama = {}; // γ

		calc_reflection_and_deltaPhase(m+1) //bad but fast :D; calc between m-1 and m
		l(' ')

		for(var k = m; k >= 2; k--){
			calc_reflection_and_deltaPhase(k);
			calc_deg_gama(k);
			calc_tg_deltaPhase_and_deltaPhase(k);
			calc_deg_beta(k);
			calc_deg_alpha(k);
			calc_multipleReflection(k);
		}

		calc_R(lambda);

	}

	//l(reflection)

	//l('Ri : ',Ri[1].toFixed(3))
	l(' ');l(' ');
	l('results : ')
	for(var key in Ri){
		l((key*1000).toFixed(0), 'нм : ', (Ri[key]*100).toFixed(3))
	}




	function calc_reflection_and_deltaPhase(k){
		var index = '' + (k-2) + (k-1),
			n1 = layers[k-2].n,
			n2 = layers[k-1].n;


		reflection[index] = Math.abs( (n1 - n2) / (n1 + n2));
		//reflection[index] = (n1 - n2) / (n1 + n2);
		deltaPhase[index] = n2 < n1 ? 0 : PI;

		l(k, ' : ', k-2, k-1)
		l(`reflection[${index}] : `, reflection[index].toFixed(3))
		l(`deltaPhase[${index}] : `, (deltaPhase[index] / Math.PI).toFixed(1) + ' pi')
	}


	function calc_deg_gama(k){
		var deltaPhaseIndex = '' + (k-1) + m; // ∆(k-1),m
		deg_gama[k - 2] = deltaPhase[deltaPhaseIndex] - 4 * PI * layers[k-1].n * layers[k-1].h / lambda;

		
		l(`deltaPhase[${deltaPhaseIndex}] : ${(deltaPhase[deltaPhaseIndex]/Math.PI).toFixed(1)} pi`)
		l(`deg_gama[${k-2}] : ${(deg_gama[k-2]/Math.PI).toFixed(1)} pi`);
	}

	function calc_tg_deltaPhase_and_deltaPhase(k){
		//indexes :
		// k2m -> k-2,m ...
		var k2m = '' + (k-2) + m,
			k1m = '' + (k-1) + m,
			k2k1 = '' + (k-2) + (k-1),
			k2 = '' + (k-2);

		//elements
		var r_k1m = reflection[k1m],
			r_k2k1 = reflection[k2k1],
			d_k2k1 = deltaPhase[k2k1],
			deg_g_k2 = deg_gama[k2],
			//dont go to degree
			//sin_g = Math.sin(DegToRad(deg_g_k2)),
			//cos_g = Math.cos(DegToRad(deg_g_k2)),
			//cos_d = Math.cos(DegToRad(d_k2k1));
			sin_g = Math.sin(deg_g_k2),
			cos_g = Math.cos(deg_g_k2),
			cos_d = Math.cos(d_k2k1);


		var tg = (r_k1m * (1 - Math.pow(r_k2k1, 2)) * sin_g) 
				/ (r_k2k1 * (1 + Math.pow(r_k1m, 2)) *cos_d + r_k1m*(1 + Math.pow(r_k2k1, 2)) * cos_g);

		//dont go to degree
		//deltaPhase[k2m] = RadToDeg(Math.atan(tg));
		deltaPhase[k2m] = Math.atan(tg);


		l(`tg : ${tg.toFixed(3)}`);
		l(`arctg : ${(Math.atan(tg) / Math.PI).toFixed(3)} pi`)

		//l(`deltaPhase[${k2m}] : ${deltaPhase[k2m]}`)		
	}


	function calc_deg_beta(k){
		//indexes :
		// k2m -> k-2,m ...
		var k1 = '' + (k-1),
			k2 = '' + (k-2),
			k1m = '' + (k-1) + m,
			k2k1 = '' + (k-2) + (k-1),
			n = layers[k1].n,
			h = layers[k1].h;
		
		var beta = deltaPhase[k2k1] + deltaPhase[k1m] - (4 * PI * n * h)/lambda;

		deg_beta[k2] = beta;

		l(`beta : ${(beta / Math.PI).toFixed(1)} pi`)
	}


	function calc_deg_alpha(k){
		//indexes :
		// k2m -> k-2,m ...
		var k1 = '' + (k-1),
			k2 = '' + (k-2),
			k1m = '' + (k-1) + m,
			k2k1 = '' + (k-2) + (k-1),
			n = layers[k1].n,
			h = layers[k1].h;
		
		var alpha = -deltaPhase[k2k1] + deltaPhase[k1m] - (4 * PI * n * h)/lambda;

		deg_alpha[k2] = alpha;

		l(`alpha : ${(alpha / Math.PI).toFixed(1)} pi`)
	}


	function calc_multipleReflection(k){
		//indexes :
		// k2m -> k-2,m ...
		var k2 = '' + (k-2),
			k2m = '' + (k-2) + m,
			k1m = '' + (k-1) + m,
			k2k1 = '' + (k-2) + (k-1);

		//elements
		var	r_k2k1 = reflection[k2k1],
			r2_k2k1 = Math.pow(r_k2k1, 2),
			r_k1m = reflection[k1m],
			r2_k1m = Math.pow(r_k1m, 2),
			//dont go to degree
			//cos_a = Math.cos(DegToRad( deg_alpha[k2] )),
			//cos_b = Math.cos(DegToRad( deg_beta[k2] )),
			cos_a = Math.cos(deg_alpha[k2]),
			cos_b = Math.cos(deg_beta[k2]),
			reflection2;

		reflection2 = (r2_k2k1 + r2_k1m + 2 * r_k2k1 * r_k1m * cos_a) /
						  (1 + r2_k2k1 * r2_k1m + 2 * r_k2k1 * r_k1m * cos_b);
		reflection[k2m] = Math.sqrt(reflection2)						  

		//l('deg_alpha[k2] : ',(deg_alpha[k2] / Math.PI).toFixed(3) + ' pi')
		//l('deg_beta[k2] : ',(deg_beta[k2] / Math.PI).toFixed(3) + ' pi')
		l('cos_a : ', cos_a)
		l('cos_b : ', cos_b)
		
		l(`r_k1m(r${k1m}) : `, r_k1m.toFixed(3))
		l(`r2_k1m(r2${k1m}) : `, r2_k1m.toFixed(3))
		l('r_k2k1 : ', r_k2k1.toFixed(3))
		l('r2_k2k1 : ', r2_k2k1.toFixed(3))

		l(`reflection[${k2m}] : `, reflection[k2m].toFixed(3))	
		l(' ')			  
	}


	function calc_R(lambda){
		//index
		var index = '0' + m 

		Ri[lambda] = reflection[index] * reflection[index];
	}



	//extra functions
	function DegToRad(degree){
		return degree * Math.PI / 180
	}

	function RadToDeg(rad){
		return rad * 180 / Math.PI
	}
}