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
  			<h4 className="description-header">Description</h4>
  			<p>RoLiM iteratively detects over-represented linear motifs in sequence data sets.</p>
  		</section>
  		<JobForm endpoint="/api/job/" />
  	</main>
  	<footer>
  		<p>Copyright &copy; 2020 Theodore G. Smith</p>
  	</footer>
  </div>
);
const wrapper = document.getElementById("app");
wrapper ? ReactDOM.render(<App />, wrapper) : null;
