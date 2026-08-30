# *BAGS: Barcode, Audit & Grade System*

<div style="border: 2px solid red; padding: 15px; background-color: #ffe6e6; font-size: 18px;">
  <h1 style="color:red; text-align:center;">⚠️ Important Notice ⚠️</h1>

  <p>
    <strong>BAGS is currently not operational.</strong>
  </p>

  <p>
    BAGS was originally developed to retrieve barcode records from BOLD Systems using the BOLD v4 API and the R package <code>bold</code>. Following the transition of BOLD Systems to version 5 and its new API infrastructure, the v4 services on which BAGS relies no longer provide the data access required by the application. Consequently, BAGS is currently unable to retrieve the BOLD records needed to perform its analyses and grade assignments.
  </p>

  <p>
    In addition, the R package <code>bold</code>, which BAGS uses for data retrieval, was removed from the CRAN repository in August 2024. Although archived versions of the package remain available, they were developed for the BOLD v4 API and therefore do not restore the current functionality of BAGS.
  </p>

 <p> In light of this, we recommend <strong>BOLDcuratoR</strong>, a tool developed within the <strong>Biodiversity Genomics Europe (BGE)</strong> project for curating barcode data retrieved through the current BOLD infrastructure. Among its features, BOLDcuratoR includes the BAGS A–E grade assignment, alongside additional tools for assessing BIN concordance, ranking specimens, and supporting reference-library curation. </p>

<p> The BOLDcuratoR web application is available at <a href="https://benprice.shinyapps.io/BOLDcuratoR/" target="_blank">https://benprice.shinyapps.io/BOLDcuratoR/</a>, and its source code is available on <a href="https://github.com/bge-barcoding/BOLDcuratoR" target="_blank">GitHub</a>. </p>
