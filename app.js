const buttonCalcul = document.getElementById("button-calcul");
// future imputs
const q_c_num = document.getElementById("q_c");
const w_p_num = document.getElementById("w_p");
const b_num = document.getElementById("b");
const debit_agent_num = document.getElementById("debit_agent");
const rho_apa_num = document.getElementById("rho_apa");
const rho_agentG_num = document.getElementById("rho_agentG");
const rho_agentL_num = document.getElementById("rho_agentL");
const miu_apa_num = document.getElementById("miu_apa");
const miu_agentG_num = document.getElementById("miu_agentG");
const lambda_apa_num = document.getElementById("lambda_apa");
const lambda_num = document.getElementById("lambda");
const lambda_agent_num = document.getElementById("lambda_agent");
const lambda_placa_num = document.getElementById("lambda_placa");
const t_c_num = document.getElementById("t_c");
const t_wi_num = document.getElementById("t_wi");
const t_we_num = document.getElementById("t_we");
const miu_p_num = document.getElementById("miu_p");
const h_i_num = document.getElementById("h_i");
const h_e_num = document.getElementById("h_e");
const l_p_num = document.getElementById("l_p");
const delta_p_num = document.getElementById("delta_p");
let beta_num = document.getElementById("beta");

//rezultate
const n_p_span = document.getElementById("n_p");
const deltap_apa_span = document.getElementById("deltap_apa");
const deltap_agent_span = document.getElementById("deltap_agent");
const w_agentG_span = document.getElementById("w_agentG");
const w_agentL_span = document.getElementById("w_agentL");
const d_i_agent_span = document.getElementById("d_i_agent");
const d_e_agent_span = document.getElementById("d_e_agent");
const d_i_apa_span = document.getElementById("d_i_apa");
const w_span = document.getElementById("w");
const eroare_span = document.getElementById("eroare");
const iteratii_span = document.getElementById("iteratii");

// Valori de intrare si constante
let q_c;
let w_p;
let b;
let debit_agent;
let rho_apa;
let rho_agentG;
let rho_agentL;
let miu_apa;
let miu_agentG;
let lambda_apa;
let lambda;
let lambda_agent;
let lambda_placa;
let t_c;
let t_wi;
let t_we;
let miu_p;
let h_i;
let h_e;
let l_p;
let delta_p;
let beta;

let a;
let n_p_nerot;
let n_p;
let deltap_apa;
let a_c;
let n_c;
let n_c_rot;
let w_agent;
let w_apa;
let re_apa;
let pr_apa;
let p;
let f_1;
let f_t;
let f_a;
let m;
let nu_1;
let nu_t;
let nu_a;
let alfa_apa;
let alfa_agent;
let re_agent;
let f_f;
let deltap_agent;
let w_agentG;
let w_agentL;
let d_i_agent;
let d_e_agent;
let d_i_apa;
let d_e_apa;
let w;
let a_p;
let x;
let theta;
let d_h;
let c_pv;
let r_f;
let r_a;
let debit_apa;
let deltat_m;
let f;
let deltat_m_inf;
let t_med;
let t_p;
let g;
let k_0;
let k_calc;
let eroare;
let iteratii;
// Calcule necesare
function calcule() {

  a_p = l_p * w_p;
  x = b * Math.PI / lambda;
  theta = 1 / 6 * (1 + Math.sqrt(1 + x ** 2) + 4 * Math.sqrt(1 + x ** 2 / 2));
  d_h = 2 * b / theta;
  c_pv = 4.175;
  r_f = 1.75 * 10 ** (-4);
  r_a = 2.17 * 10 ** (-4);
  debit_apa = q_c / (c_pv * (t_we - t_wi));
  deltat_m = ((t_c - t_we) - (t_c - t_wi)) / Math.log((t_c - t_we) / (t_c - t_wi));
  f = 0.967;
  deltat_m_inf = f * deltat_m;
  t_med = (30 + 34) / 2;

  t_p = (t_med + t_c) / 2;

  g = 9.81;
  k_0 = 500;
  k_calc = 100;
  eroare = 80;
  iteratii = 0;

  while (eroare > 5) {
    // Suprafata de caldura
    a = (q_c * 1000) / (k_0 * deltat_m);

    // Numărul total de plăci
    n_p_nerot = a / a_p;
    n_p = Math.ceil(n_p_nerot);

    // Aria transversală a canalului de curgere
    a_c = w_p * b;

    // Numărul de canale de trecere pentru fiecare fluid
    n_c = (n_p - 1) / 2;
    n_c_rot = Math.ceil(n_c);
    if (n_c_rot % 2 == 1) {
      f = 0.42;
    } else {
      f = 0.967;
    }

    // Astfel. diferenţa de temperatură medie logaritmică corectată se recalculează
    deltat_m_inf = f * deltat_m;

    // Viteza de curgere a agentului frigorific prin condensator
    w_agent = debit_agent / (rho_agentL * a_c * n_c_rot);

    // Viteza de curgere a apei de răcire prin condensator
    w_apa = debit_apa / (rho_apa * a_c * n_c_rot);

    // Numărul Reynolds pentru apă
    re_apa = (w_apa * d_h * rho_apa) / miu_apa;

    // Numărul Prandtl pentru apă
    pr_apa = c_pv * miu_apa / lambda_apa;

    // Determinarea factorului de frecare Fanning

    p = 0.00423 * beta + 0.0000223 * beta ** 2;

    // Componenta longitudinală a factorului de frecare Fanning
    f_1 = 1774 * beta ** (-1.026) * theta ** (2) * re_apa ** (-1);

    // Componenta transversală a factorului de frecare Fanning
    f_t = 46.6 * beta ** (-1.08) * theta ** (1 + p) * re_apa ** (-p);

    // Factorul de frecare Fanning
    f_a = (f_1 ** 3 + f_t ** 3) ** (1 / 3);

    // Numărul  Nusselt se calculează pe direcţie longitudinală şi transversală
    m = 0.646 + 0.0011 * beta;

    // Numărul Nusselt pe direcţie longitudinală
    nu_1 = 3.65 * beta ** (-0.455) * theta ** (0.661) * re_apa ** (0.339);

    //Numărul Nusselt pe direcţie transversală
    nu_t = 12.6 * beta ** (-1.142) * theta ** (1 - m) * re_apa ** (m);

    //Numărul Nusselt pentru apă
    nu_a = (nu_1 ** 3 + nu_t ** 3) ** (1 / 3) * pr_apa ** (1 / 3) * (miu_apa / miu_p) ** (0.17);

    // Coeficientul de transfer de căldură prin convecţie între placă şi apă
    alfa_apa = (lambda_apa * nu_a) / d_h;

    // Cunoscând entalpia agentului frigorific la intrarea şi la ieşirea din condensator.
    // se poate calcula căldură latent de condensare agentul frigorific
    r = h_i - h_e;

    // Coeficientul de transfer de căldură prin convecţie între placă şi agentul frigorific
    alfa_agent = 0.943 * (((rho_agentL * (rho_agentG + rho_agentL) * g * r * 1000 * lambda_agent) /
      (miu_agentG * l_p * (t_c - t_p))) ** (0.25));

    // Numărul Reynold pentru agentul frigorific
    re_agent = (w_agent * d_h * rho_agentG) / miu_agentG;

    // Factorul de frecare
    f_f = 0.6 * (re_agent ** (-0.3));

    // Coeficientul global de transfer de căldură
    k_calc = 1 / (1 / alfa_agent + 1 / alfa_apa + r_f + r_a + delta_p / lambda_placa);

    // Eroarea de calcul
    eroare = Math.abs((k_calc - k_0) / k_0) * 100;

    k_0 = k_calc;
    k_calc = 100;

    // numar iteratii
    iteratii = iteratii + 1;

  }

  // Pierderile de presiune din condensator
  deltap_apa = f_a * rho_apa * w_apa ** 2 * 2 * (l_p / d_h);
  deltap_agent = 8 * f_f * rho_agentG * w_agent ** 2 / 2 * (l_p / d_h);

  // Dimensiunile geometrice ale condensatorului
  //Viteza de curgere a agentului frigorific prin condensator
  //la intrare:
  w_agentG = debit_agent / (rho_agentG * a_c * n_c);
  // la iesire
  w_agentL = debit_agent / (rho_agentL * a_c * n_c);

  // Diametrele conductelor pentru agentul frigorific sunt
  // la intrare
  d_i_agent = Math.sqrt(2 * debit_agent / (Math.PI * rho_agentG * w_agentG));
  // la iesire
  d_e_agent = Math.sqrt(2 * debit_agent / (Math.PI * rho_agentL * w_agentL));

  // Viteza de curgere a apei de răcire prin condensator se consideră constantă
  // Diametrul conductei de intrare a apei în condensator
  d_i_apa = Math.sqrt(2 * debit_apa / (Math.PI * rho_apa * w_apa));
  // Deoarece apa nu îşi modifică starea de agregare, iar temperatură creşte în mod
  // nesemnificativ de la intrarea şi ieşirea sa din condensator
  d_e_apa = d_i_apa;

  // Lungimea W a condensatorului depinde de numărul de plăci al condesatorului şi de
  // grosimea unei plăci:
  w = n_p_nerot * delta_p;

}


q_c_num.addEventListener('input', () => {
  q_c = parseFloat(q_c_num.value);
  localStorage.setItem("q_c", q_c);
});
w_p_num.addEventListener('input', () => {
  w_p = parseFloat(w_p_num.value);
  localStorage.setItem("w_p", w_p);
});
b_num.addEventListener('input', () => {
  b = parseFloat(b_num.value);
  localStorage.setItem("b", b);
});
debit_agent_num.addEventListener('input', () => {
  debit_agent = parseFloat(debit_agent_num.value);
  localStorage.setItem("debit_agent", debit_agent);
});
rho_apa_num.addEventListener('input', () => {
  rho_apa = parseFloat(rho_apa_num.value);
  localStorage.setItem("rho_apa", rho_apa);
});
rho_agentG_num.addEventListener('input', () => {
  rho_agentG = parseFloat(rho_agentG_num.value);
  localStorage.setItem("rho_agentG", rho_agentG);
});
rho_agentL_num.addEventListener('input', () => {
  rho_agentL = parseFloat(rho_agentL_num.value);
  localStorage.setItem("rho_agentL", rho_agentL);
});
miu_apa_num.addEventListener('input', () => {
  miu_apa = parseFloat(miu_apa_num.value);
  localStorage.setItem("miu_apa", miu_apa);
});
miu_agentG_num.addEventListener('input', () => {
  miu_agentG = parseFloat(miu_agentG_num.value);
  localStorage.setItem("miu_agentG", miu_agentG);
});
lambda_apa_num.addEventListener('input', () => {
  lambda_apa = parseFloat(lambda_apa_num.value);
  localStorage.setItem("lambda_apa", lambda_apa);
});
lambda_num.addEventListener('input', () => {
  lambda = parseFloat(lambda_num.value);
  localStorage.setItem("lambda", lambda);
});
lambda_agent_num.addEventListener('input', () => {
  lambda_agent = parseFloat(lambda_agent_num.value);
  localStorage.setItem("lambda_agent", lambda_agent);
});
lambda_placa_num.addEventListener('input', () => {
  lambda_placa = parseFloat(lambda_placa_num.value);
  localStorage.setItem("lambda_placa", lambda_placa);
});
t_c_num.addEventListener('input', () => {
  t_c = parseFloat(t_c_num.value);
  localStorage.setItem("t_c", t_c);
});
t_wi_num.addEventListener('input', () => {
  t_wi = parseFloat(t_wi_num.value);
  localStorage.setItem("t_wi", t_wi);
});
t_we_num.addEventListener('input', () => {
  t_we = parseFloat(t_we_num.value);
  localStorage.setItem("t_we", t_we);
});
miu_p_num.addEventListener('input', () => {
  miu_p = parseFloat(miu_p_num.value);
  localStorage.setItem("miu_p", miu_p);
});
h_i_num.addEventListener('input', () => {
  h_i = parseFloat(h_i_num.value);
  localStorage.setItem("h_i", h_i);
});
h_e_num.addEventListener('input', () => {
  h_e = parseFloat(h_e_num.value);
  localStorage.setItem("h_e", h_e);
});
l_p_num.addEventListener('input', () => {
  l_p = parseFloat(l_p_num.value);
  localStorage.setItem("l_p", l_p);
});
delta_p_num.addEventListener('input', () => {
  delta_p = parseFloat(delta_p_num.value);
  localStorage.setItem("delta_p", delta_p);
});
beta_num.addEventListener('input', () => {
  beta = parseFloat(beta_num.value);
  localStorage.setItem("beta", beta);
});

//local storage
q_c_num.value = localStorage.getItem('q_c');
w_p_num.value = localStorage.getItem("w_p");
b_num.value = localStorage.getItem("b");
debit_agent_num.value = localStorage.getItem("debit_agent");
rho_apa_num.value = localStorage.getItem("rho_apa");
rho_agentG_num.value = localStorage.getItem("rho_agentG");
rho_agentL_num.value = localStorage.getItem("rho_agentL");
miu_apa_num.value = localStorage.getItem("miu_apa");
miu_agentG_num.value = localStorage.getItem("miu_agentG");
lambda_apa_num.value = localStorage.getItem("lambda_apa");
lambda_num.value = localStorage.getItem("lambda");
lambda_agent_num.value = localStorage.getItem("lambda_agent");
lambda_placa_num.value = localStorage.getItem("lambda_placa");
t_c_num.value = localStorage.getItem("t_c");
t_wi_num.value = localStorage.getItem("t_wi");
t_we_num.value = localStorage.getItem("t_we");
miu_p_num.value = localStorage.getItem("miu_p");
h_i_num.value = localStorage.getItem("h_i");
h_e_num.value = localStorage.getItem("h_e");
l_p_num.value = localStorage.getItem("l_p");
delta_p_num.value = localStorage.getItem("delta_p");
beta_num.value = localStorage.getItem("beta");

q_c = parseFloat(q_c_num.value);
w_p = parseFloat(w_p_num.value);
b = parseFloat(b_num.value);
debit_agent = parseFloat(debit_agent_num.value);
rho_apa = parseFloat(rho_apa_num.value);
rho_agentG = parseFloat(rho_agentG_num.value);
rho_agentL = parseFloat(rho_agentL_num.value);
miu_apa = parseFloat(miu_apa_num.value);
miu_agentG = parseFloat(miu_agentG_num.value);
lambda_apa = parseFloat(lambda_apa_num.value);
lambda = parseFloat(lambda_num.value);
lambda_agent = parseFloat(lambda_agent_num.value);
lambda_placa = parseFloat(lambda_placa_num.value);
t_c = parseFloat(t_c_num.value);
t_wi = parseFloat(t_wi_num.value);
t_we = parseFloat(t_we_num.value);
miu_p = parseFloat(miu_p_num.value);
h_i = parseFloat(h_i_num.value);
h_e = parseFloat(h_e_num.value);
l_p = parseFloat(l_p_num.value);
delta_p = parseFloat(delta_p_num.value);
beta = parseFloat(beta_num.value);

//functia de calcul
function valori() {
  calcule()
  n_p_span.innerHTML = n_p;
  deltap_apa_span.innerHTML = deltap_apa.toFixed(5);
  deltap_agent_span.innerHTML = deltap_agent.toFixed(10);
  w_agentG_span.innerHTML = w_agentG.toFixed(4);
  w_agentL_span.innerHTML = w_agentL.toFixed(10);
  d_i_agent_span.innerHTML = d_i_agent.toFixed(10);
  d_e_agent_span.innerHTML = d_e_agent.toFixed(10);
  d_i_apa_span.innerHTML = d_i_apa.toFixed(10);
  w_span.innerHTML = w.toFixed(5);
  eroare_span.innerHTML = eroare;
  iteratii_span.innerHTML = iteratii;
  console.log(n_c_rot);
}

buttonCalcul.addEventListener("click", valori);
