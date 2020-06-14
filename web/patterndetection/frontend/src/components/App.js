import React from "react";
import ReactDOM from "react-dom";
import DataProvider from "./DataProvider";
import Table from "./Table";
import JobForm from "./JobForm";

const App = () => (
  <div className="outer">
  	<header>
  		<div className="labtitle">
  			<div className="labtitle-inner">
  				<h1><a className="lablink" href="https://langelab.med.ubc.ca/">The Lange Lab - for Translational Proteomics in Childhood Cancer</a></h1>
  			</div>
  		</div>
  		<div className="tooltitle">
  			<div className="tooltitle-inner">
  				<h2 className="title">RoLiM</h2>
  				<h3 className="subtitle">Robust detection of linear motifs in sequence data.</h3>
  			</div>
  		</div>
  	</header>
  	<br />
  	<main>
  		<section>
  			<h4 className="title is-4 description-header">Description</h4>
  			<p className="description">RoLiM iteratively detects over-represented linear motifs in sequence data sets.</p>
  		</section>
  		<section>
  			<h4 className="title is-4 description-header">Citing RoLiM</h4>
  			<p className="description">Please visit our <a href="https://langelab.med.ubc.ca">main lab website</a> for more information about citing RoLiM.</p>
  		</section>
  		<JobForm endpoint="/api/job/" />
  	</main>
  	<footer>
  		<div className="footer-inner">
  			<p className="copyright">Copyright &copy; 2020 Theodore G. Smith</p>
  		</div>
  	</footer>
  </div>
);
const wrapper = document.getElementById("app");
wrapper ? ReactDOM.render(<App />, wrapper) : null;
